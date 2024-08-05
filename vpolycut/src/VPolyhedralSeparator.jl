using SCIP
using JuMP
using LinearAlgebra
import MathOptInterface as MOI

"""
Intersection Cut Separator
# Fields
- scipd::SCIP.SCIPData Reference to the SCIPData object
"""
@kwdef mutable struct VPolyhedralSeparator <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
    called::Int = 0
end

function get_point_ray_collection(
    scip::SCIP.SCIPData, root_tableau::Tableau, path::Vector{Node}
)::Tuple{Vector{Point},Vector{Ray},Vector{SCIP.SCIP_Real}}
    # Enter Probing mode
    SCIP.SCIPstartProbing(scip)

    # Get To Node
    for node in path
        if isroot(node)
            continue
        end
        do_action(scip, get_action(node))
    end

    # Propegate
    prunable = propagate!(scip)
    if prunable
        SCIP.SCIPendProbing(scip)
        return [], [], []
    end

    # Solve LP
    lp_feasible = solve_lp_relaxation(scip)
    if !lp_feasible
        SCIP.SCIPendProbing(scip)
        return [], [], []
    end

    # Get Optimal Tableau
    tableau = construct_tableau(scip)
    corner = construct_corner_polyhedron(tableau)
    solution = SCIP.SCIPgetLPObjval(scip)
    @assert get_nvars(tableau) == get_nvars(root_tableau)
    # Leave Probing mode
    SCIP.SCIPendProbing(scip)
    return [get_lp_sol(corner)], get_lp_rays(corner), [solution]
end

function solve_separating_lp(lp_solution, intersection_points, pararrel_rays)
    dim = length(lp_solution)
    separating_lp = create_subscip_model()

    @assert SCIP.SCIPgetSubscipsOff(unsafe_backend(separating_lp).inner) != 0

    @variable(separating_lp, x[1:dim])

    for point in intersection_points
        new_point = point - lp_solution
        @constraint(separating_lp, sum(x[i] * new_point[i] for i in 1:dim) >= 1)
    end

    for ray in pararrel_rays
        @constraint(separating_lp, sum(x[i] * ray[i] for i in 1:dim) >= 0)
    end

    @objective(separating_lp, Min, 0)
    try
        optimize!(separating_lp)
    catch e
        str = randstring(5)
        scip = unsafe_backend(separating_lp).inner
        println("Error in solving seperating LP written $(str)")
        SCIP.@SCIP_CALL SCIP.SCIPwriteLP(scip, "separating_lp_$(str).lp")
        return nothing
    end
    write_to_file(separating_lp, "separating_lp.lp")
    display(solution_summary(separating_lp))
    if is_solved_and_feasible(separating_lp)
        return Point(value.(x))
    end

    return nothing
end

"""
Construct the seperating lp return a reference to a pointer to the LPI
"""
function find_cut_from_split(
    sepa::VPolyhedralSeparator,
    split_index::Int64,
    root_corner_polyhedron::CornerPolyhedron,
    tableau::Tableau
)::Union{Nothing,Ref{Ptr{SCIP.SCIP_ROW}}}
    scip = sepa.scipd
    lp_sol = corner_polyhedron.lp_sol
    lp_rays = corner_polyhedron.lp_rays

    # Setup Point Ray Gatherer
    println("Starting Branch and Bound")
    split_var = get_var_from_column(tableau, split_index)
    split_var_ptr = get_var_pointer(split_var)
    brachandbound = BranchAndBound(
        scip; max_leaves=2, branching_rule=PriorityBranching(split_var_ptr)
    )
    execute_branchandbound(brachandbound)

    # STEP 1: Get Point Ray
    intersection_points = []
    parallel_ray = []
    s = 1
    lp_sol = get_lp_sol(root_corner_polyhedron)
    for leaf in get_leaves(brachandbound)
        s += 1
        path = get_path(leaf)
        points, rays = get_point_ray_collection(scip, tableau, path)
        append!(intersection_points, points)
        append!(parallel_ray, rays)
    end
    println(
        "Point ray collection contains $(length(intersection_points)) points and $(length(parallel_ray)) rays"
    )

    # Create a projection object 
    projection = create_projection_to_nonbasic_space(tableau)

    # STEP 2: Project the intersection points and rays to the non-basic space
    projected_lp_sol = project_point(projection, lp_sol)
    projected_intersection_points = [
        project_point(projection, point) for point in intersection_points
    ]
    projected_parallel_ray = [project_point(projection, ray) for ray in parallel_ray]
    @info "Projected Sol" [projected_lp_sol]
    @info "Projected Intersection" projected_intersection_points
    @info "Projected Ray" projected_parallel_ray

    # Step 3: Solve the seperating LP
    separating_sol = solve_separating_lp(
        projected_lp_sol, projected_intersection_points, projected_parallel_ray
    )
    @info "Seperating Solution" [separating_sol]
    if isnothing(separating_sol)
        return nothing
    end

    # Step 4: Convert the seperating solution to the original space
    full_seperating_sol = undo_projection(projection, separating_sol)
    # Undo Complementation
    full_uncomplemented_sol = get_uncomplemented_vector(
        full_seperating_sol, complemented_tableau
    )
    uncomplmented_lp_sol = get_uncomplemented_vector(lp_sol, complemented_tableau)
    println(full_uncomplemented_sol)
    b = dot(full_uncomplemented_sol, uncomplmented_lp_sol) + 1

    # Step 5: Convert the seperating solution to a cut in general form
    # This uses tableau rows from constraint_matrix so complementaztion does not matter
    cut_vector, b = convert_standard_inequality_to_general(
        scip, tableau, full_uncomplemented_sol, b
    )

    # Step 6: Add the cut to the SCIP
    row = Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRowSepa(
        scip,
        row,
        scip.sepas[sepa],
        "",
        b,
        SCIP.SCIPinfinity(scip),
        true,
        false,
        false
    )
    infeasible = Ref{SCIP.SCIP_Bool}(0)

    for (idx, sol) in enumerate(cut_vector)
        if SCIP.SCIPisZero(scip, sol) == 0
            var = get_var_from_column(tableau, idx)
            SCIP.@SCIP_CALL SCIP.SCIPaddVarToRow(scip, row[], get_var_pointer(var), sol)
        end
    end

    #SCIP.@SCIP_CALL SCIP.SCIPprintRow(scip, row[], C_NULL)
    return row
    SCIP.@SCIP_CALL SCIP.SCIPaddRow(scip, row[], true, infeasible)
    SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip, row)
    return true
end

function SCIP.exec_lp(sepa::VPolyhedralSeparator)
    # Aliasing for easier call
    scip = sepa.scipd
    sepa.called += 1

    # Preconditions
    @assert(SCIP.SCIPgetStage(scip) == SCIP.SCIP_STAGE_SOLVING)
    @assert(SCIP.SCIPisLPSolBasic(scip) != 0)
    @assert(SCIP.SCIPgetLPSolstat(scip) == SCIP.SCIP_LPSOLSTAT_OPTIMAL)

    # STEP 0: Get LP Tableau data
    tableau = construct_tableau_with_constraint_matrix(scip)

    # STEP 1: Get Corner Polyhedron of current LP solution
    corner = construct_corner_polyhedron(tableau)
    @info "Corner Polyhedron" [corner.lp_sol] corner.lp_rays
    # STEP 2: Get the set of fractional indices 
    split_indices = get_branching_indices(scip, tableau)
    pool::Vector{Ref{Ptr{SCIP.SCIP_ROW}}} = []

    # STEP 3: For each fractional indices try to find a cut
    for (i, index) in enumerate(split_indices)
        cut = find_cut_from_split(sepa, index, corner, tableau)
        if !isnothing(cut)
            push!(pool, cut)
        end
    end
    println("Found $(length(pool)) cuts")
    for cut in pool
        infeasible = Ref{SCIP.SCIP_Bool}(0)
        SCIP.SCIPprintRow(scip, cut[], C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPaddRow(scip, cut[], true, infeasible)
        SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip, cut)
    end

    if length(pool) > 0
        return SCIP.SCIP_SEPARATED
    end

    return SCIP.SCIP_DIDNOTFIND
end

function include_vpolyhedral_sepa(scip::SCIP.SCIPData)
    sepa = VPolyhedralSeparator(; scipd=scip)
    SCIP.include_sepa(scip.scip[], scip.sepas, sepa; freq=0, usessubscip=true)
end
