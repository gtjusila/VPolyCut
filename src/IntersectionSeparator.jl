using SCIP
using JuMP
using LinearAlgebra
import MathOptInterface as MOI

"""
Intersection Cut Separator
# Fields
- scipd::SCIP.SCIPData Reference to the SCIPData object
"""
@kwdef mutable struct IntersectionSeparator <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
    called::Int = 0
end

"""
Given a current solution, a split index, and a ray collection compute intersection points
and pararrel_rays
"""
function compute_intersection_points(
    scip::SCIP.SCIPData,
    split_index::Int,
    current_solution::Point,
    ray_collection::Vector{Ray}
)
    intersection_points = []
    pararrel_ray = []

    for ray in ray_collection
        low = floor(current_solution[split_index])
        up = ceil(current_solution[split_index])
        epsilon = 0

        if is_negative(scip, ray[split_index])
            # Ray will hit lower bound
            epsilon = (low - current_solution[split_index]) / ray[split_index]
        elseif is_positive(scip, ray[split_index])
            # Ray will hit upper bound
            epsilon = (up - current_solution[split_index]) / ray[split_index]
        else
            # Ray is parallel_ray to the split
            push!(pararrel_ray, ray)
            continue
        end

        # Ray will intersect bound so add to intersection point set
        push!(intersection_points, current_solution + epsilon * ray)
    end

    return (intersection_points, pararrel_ray)
end

function solve_intersection_separating_lp(lp_solution, intersection_points, pararrel_rays)
    dim = length(lp_solution)

    separating_lp = setup_scip_safe_jump_model()
    set_attribute(separating_lp, "display/verblevel", 0)
    set_attribute(separating_lp, "presolving/maxrounds", 0)
    scip = get_scip_data_from_model(separating_lp)
    SCIP.@SCIP_CALL SCIP.SCIPsetSubscipsOff(scip, true)

    @assert SCIP.SCIPgetSubscipsOff(unsafe_backend(separating_lp).inner) != 0

    @variable(separating_lp, x[1:dim])
    for point in intersection_points
        new_point = point - lp_solution
        @constraint(separating_lp, sum(x[i] * new_point[i] for i in 1:dim) == 1)
    end

    for ray in pararrel_rays
        @constraint(separating_lp, sum(x[i] * ray[i] for i in 1:dim) == 0)
    end

    @objective(separating_lp, Min, 0)

    optimize!(separating_lp)
    println(separating_lp)
    if is_solved_and_feasible(separating_lp)
        return Point(value.(x))
    end

    return nothing
end

"""
Construct the seperating lp return a reference to a pointer to the LPI
"""
function find_cut_from_split(
    sepa::IntersectionSeparator,
    split_index::Int64,
    corner_polyhedron::CornerPolyhedron,
    tableau::Tableau
)::Bool
    scip = sepa.scipd
    complemented_tableau = copy_and_complement(tableau)
    corner = construct_corner_polyhedron(complemented_tableau)
    lp_sol = corner.lp_sol
    lp_rays = corner.lp_rays
    tableau = get_tableau(complemented_tableau)
    @info "Complemented Corner Polyhedron" [lp_sol] lp_rays complemented_tableau.complemented_columns

    # STEP 1: Compute intersection points
    intersection_points, parallel_ray = compute_intersection_points(
        scip, split_index, lp_sol, lp_rays
    )
    parallel_ray = lp_rays
    @info "Intersection Points:" intersection_points
    @info "Rays" parallel_ray

    # Create a projection object 
    projection = create_projection_to_nonbasic_space(tableau)

    # STEP 2: Project the intersection points and rays to the non-basic space
    projected_lp_sol = project(projection, lp_sol)
    projected_intersection_points = [
        project(projection, point) for point in intersection_points
    ]
    projected_parallel_ray = [project(projection, ray) for ray in parallel_ray]
    @info "Projected Sol" [projected_lp_sol]
    @info "Projected Intersection" projected_intersection_points
    @info "Projected Ray" projected_parallel_ray

    # Step 3: Solve the seperating LP
    separating_sol = solve_intersection_separating_lp(
        projected_lp_sol, projected_intersection_points, projected_parallel_ray
    )
    @info "Seperating Solution" [separating_sol]
    if isnothing(separating_sol)
        return false
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
        if is_non_zero(scip, sol)
            var = get_var_from_column(tableau, idx)
            SCIP.@SCIP_CALL SCIP.SCIPaddVarToRow(scip, row[], get_var_pointer(var), sol)
        end
    end

    SCIP.@SCIP_CALL SCIP.SCIPprintRow(scip, row[], C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPaddRow(scip, row[], true, infeasible)
    SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip, row)
    return true
end

function SCIP.exec_lp(sepa::IntersectionSeparator)
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
    separated = false

    # STEP 3: For each fractional indices try to find a cut
    for (i, index) in enumerate(split_indices)
        @info "Splitting at index $(index)"
        if i % 10 == 0
            println("Cut Generated $(i)")
        end
        success = find_cut_from_split(sepa, index, corner, tableau)
        separated = separated || success
    end

    if (separated)
        return SCIP.SCIP_SEPARATED
    end

    return SCIP.SCIP_DIDNOTFIND
end

function include_intersection_sepa(scip::SCIP.SCIPData)
    sepa = IntersectionSeparator(; scipd=scip)
    SCIP.include_sepa(scip.scip[], scip.sepas, sepa; freq=0, usessubscip=true)
end