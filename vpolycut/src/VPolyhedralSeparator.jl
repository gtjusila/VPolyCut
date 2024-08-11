using SCIP
using JuMP
using LinearAlgebra
import MathOptInterface as MOI
using HiGHS: HiGHS

"""
VPolyhedral Cut Separator

Implementation based on Algorithm 1 of:
Balas, Egon, and Aleksandr M. Kazachkov. "V-polyhedral disjunctive cuts." arXiv preprint arXiv:2207.13619 (2022).
Disjunction are obtained from partial branch and bound trees

# Constructor Parameters
- scipd::SCIP.SCIPData Reference to the SCIPData object
- n_leaves::Int Number of leaves in the disjunction

and datastructures for the separator
"""
@kwdef mutable struct VPolyhedralSeparator <: SCIP.AbstractSeparator
    "Pointer to SCIP"
    scipd::SCIP.SCIPData
    "Number of leaves in the disjunction"
    n_leaves::Int = 2
    "Number of times the separator is called"
    called::Int = 0
    "Cut Limit"
    cut_limit::Int = -1
    "Is LP solution separated?"
    separated::Bool = false

    "Complemented Tableau"
    complemented_tableau::Union{Nothing,ComplementedTableau} = nothing
    "Disjunction"
    disjunction::Vector{Node} = []
    "PointRayCollection"
    point_ray_collection::Union{Nothing,PointRayCollection} = nothing
    "Projection"
    projection::Union{Nothing,Projection} = nothing
    "Cut Pool"
    cutpool::Union{Nothing,CutPool} = nothing
    "Separating Problem"
    separating_problem::Union{Nothing,AbstractModel} = nothing
end

# Include Helper
function include_vpolyhedral_sepa(scip::SCIP.SCIPData; n_leaves=2, cut_limit=-1)
    sepa = VPolyhedralSeparator(; cut_limit=cut_limit, scipd=scip, n_leaves=n_leaves)
    SCIP.include_sepa(
        scip.scip[], scip.sepas, sepa; priority=9999, freq=0, usessubscip=true
    )
end

function SCIP.exec_lp(sepa::VPolyhedralSeparator)
    # Aliasing for easier call
    scip = sepa.scipd
    sepa.called += 1

    # Check Preconditions and handle accordingly
    if SCIP.SCIPgetStage(scip) != SCIP.SCIP_STAGE_SOLVING
        return SCIP.SCIP_DIDNOTRUN
    end
    if SCIP.SCIPisLPSolBasic(scip) == 0
        return SCIP.SCIP_DELAYED
    end
    if SCIP.SCIPgetLPSolstat(scip) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        return SCIP.SCIP_DELAYED
    end

    # Call Separation Routine
    with_logger(NullLogger()) do
        vpolyhedralcut_separation(sepa)
    end

    if sepa.separated
        return SCIP.SCIP_SEPARATED
    else
        return SCIP.SCIP_DIDNOTFIND
    end
end

function vpolyhedralcut_separation(sepa::VPolyhedralSeparator)
    #
    # Preparation
    #
    scip = sepa.scipd
    # If cut limit is not set then set it to the number of branch candidates
    sepa.cut_limit = Int(SCIP.SCIPgetNLPBranchCands(scip))

    # Step 0: Get complemented tableau
    construct_complemented_tableau(sepa)

    # Step 1: Get Disjunction
    get_disjunction_by_branchandbound(sepa)

    #
    # Main Algorithm 
    #

    # Step 2: Setup projection used and get the points and rays from the disjunctions
    sepa.projection = create_projection_to_nonbasic_space(sepa.complemented_tableau)
    get_point_ray_collection(sepa)

    with_logger(ConsoleLogger()) do
        @info "Number of points: $(num_points(sepa.point_ray_collection))"
        @info "Number of rays: $(num_rays(sepa.point_ray_collection))"
    end

    # Step 3: Setup Cut Pool 
    sepa.cutpool = CutPool(; tableau=sepa.complemented_tableau, scip=scip)

    # Step 4: Solve Separation Problem
    solve_separation_problems(sepa)
    #projected_lp_sol = project_point(projection, lp_sol)
    #projected_points = [project_point(projection, point) for point in points_collection]
    #projected_rays = [project_point(projection, ray) for ray in rays_collection]
    #=
    separating_collection = Point[]
    with_logger(ConsoleLogger()) do
        separating_collection = get_cuts(
            scip,
            projected_lp_sol,
            projected_points,
            projected_rays,
            solution_values, cut_limit
        )
    end

    if length(separating_collection) == 0
        return nothing
    end
    for (idx, separating_sol) in enumerate(separating_collection)
        # Step 4: Convert the seperating solution back to the original space
        full_seperating_sol = undo_projection(projection, separating_sol)
        # Undo Complementation
        full_uncomplemented_sol = get_uncomplemented_vector(
            full_seperating_sol, complemented_tableau
        )
        uncomplmented_lp_sol = get_uncomplemented_vector(lp_sol, complemented_tableau)
        b = dot(full_uncomplemented_sol, uncomplmented_lp_sol) + 1
        cut_vector, b = convert_standard_inequality_to_general(
            scip, complemented_tableau, full_uncomplemented_sol, b
        )
        # Step 5: Add the cut to SCIP
        variable_pointers = get_problem_variables_pointers(complemented_tableau)
        add_sepa_row!(scip, sepa, -cut_vector, variable_pointers, -b)
    end
    =#
end

function construct_complemented_tableau(sepa::VPolyhedralSeparator)
    original_tableau = construct_tableau_with_constraint_matrix(sepa.scipd)
    sepa.complemented_tableau = ComplementedTableau(original_tableau)
    return sepa.complemented_tableau
end

function get_disjunction_by_branchandbound(sepa::VPolyhedralSeparator)
    @info "Starting Branch And Bound to generate disjunction"
    branchandbound = BranchAndBound(
        sepa.scipd; max_leaves=sepa.n_leaves
    )
    execute_branchandbound(branchandbound)
    sepa.disjunction = get_leaves(branchandbound)
    @info "Branch and bound completed"
end

function get_point_ray_collection(
    sepa::VPolyhedralSeparator
)
    scip = sepa.scipd
    disjunction = sepa.disjunction
    projection = sepa.projection

    # We only need to store projected points and ray
    sepa.point_ray_collection = PointRayCollection(scip; projection=projection)
    for term in disjunction
        path = get_path(term) # Get all actions from root to leaf
        get_disjunctive_term_information(
            sepa, path
        )
    end

    # Clean Up Duplicate Ray
    remove_duplicate_rays(sepa.point_ray_collection)
end

"""
The the point ray collection and solution value are collected from the disjunction 
The point ray collection will be given in the complemented space therefore the 
original root_tableau is required. Point and rays found are added to the sepa point ray collection directly
"""
function get_disjunctive_term_information(
    sepa::VPolyhedralSeparator,
    path::Vector{Node}
)
    # Enter Probing mode
    root_tableau = sepa.complemented_tableau
    scip = sepa.scipd

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
        return nothing, [], 0
    end

    # Solve LP
    lp_feasible = solve_lp_relaxation(scip)
    if !lp_feasible
        SCIP.SCIPendProbing(scip)
        return nothing, [], 0
    end

    # Get Optimal Tableau
    tableau = construct_tableau(scip)

    # Complement the same columns as the root tableau
    complemented_columns = get_complemented_columns(root_tableau)
    for i in complemented_columns
        var = get_var_from_column(tableau, i)
        complement_column(tableau, var)
    end

    corner = construct_corner_polyhedron(tableau)
    solution = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    point = get_lp_sol(corner)
    @assert get_nvars(tableau) == get_nvars(root_tableau)

    # Leave Probing mode
    SCIP.SCIPendProbing(scip)

    # Add rays to point ray collection
    add_point(sepa.point_ray_collection, point, solution)
    for ray in get_lp_rays(corner)
        add_ray(sepa.point_ray_collection, ray)
    end
end

function solve_separation_problems(sepa::VPolyhedralSeparator)
    scip = sepa.scipd
    lp_solution = get_solution_vector(sepa.complemented_tableau)
    lp_solution = project(sepa.projection, lp_solution)
    cut_limit = sepa.cut_limit
    dimension = length(lp_solution)

    separating_lp = Model(HiGHS.Optimizer)
    set_optimizer_attribute(separating_lp, "output_flag", false)
    points = get_points(sepa.point_ray_collection)

    translated_points = map(points) do point
        return CornerPoint(get_point(point) - lp_solution, get_objective_value(point))
    end
    rays = get_rays(sepa.point_ray_collection)

    # Construct PRLP0
    @variable(separating_lp, x[1:dimension] >= 0)
    for point in translated_points
        @constraint(separating_lp, sum(x[i] * point[i] for i in 1:dimension) >= 1)
    end
    for ray in rays
        coefficient = get_coefficients(ray)
        @constraint(separating_lp, sum(x[i] * coefficient[i] for i in 1:dimension) >= 0)
    end

    # First check if the LP is feasible by optimizing using all 0s objective
    @objective(separating_lp, Min, 0)
    optimize!(separating_lp)
    if !is_solved_and_feasible(separating_lp)
        error("Separating LP is infeasible")
    end

    # Now optimize with all 1s objective
    @objective(separating_lp, Min, sum(x))
    optimize!(separating_lp)
    if is_solved_and_feasible(separating_lp)
        cut = get_cut_from_separating_solution(sepa, value.(x))
        push!(sepa.cutpool, cut)
    else
        error("All ones objective is unbounded")
    end

    # Now Optimize with p as objective
    p_star = argmin(point -> get_objective_value(point), translated_points)
    @objective(separating_lp, Min, sum(x[i] * p_star[i] for i in 1:dimension))
    optimize!(separating_lp)
    write_to_file(separating_lp, "separating_lp_p*2.lp")
    if is_solved_and_feasible(separating_lp)
        cut = get_cut_from_separating_solution(sepa, value.(x))
        push!(sepa.cutpool, cut)
    else
        error("p* as objective is unbounded")
    end

    if !is_EQ(scip, objective_value(separating_lp), 1.0)
        @error "Cannot find cut tight at p_star is not 1.0"
    end

    # Now transition to PRLP= 
    a_bar = value.(x)
    @constraint(separating_lp, sum(x[i] * p_star[i] for i in 1:dimension) == 1)

    r_bar = filter(ray -> !is_zero(scip, dot(a_bar, ray)), rays)
    sort!(r_bar; by=ray -> abs(get_obj(get_generating_variable(ray))))

    marked = fill(false, length(r_bar))
    for ray in r_bar
        @objective(
            separating_lp, Min, sum(x[i] * ray[i] for i in 1:dimension)
        )
        optimize!(separating_lp)
        if is_solved_and_feasible(separating_lp)
            cut = get_cut_from_separating_solution(sepa, value.(x))
            push!(sepa.cutpool, cut)
        else
            with_logger(ConsoleLogger()) do
                @error "Failed to find a cut"
            end
        end

        if length(sepa.cutpool) >= cut_limit
            break
        end
    end
    add_all_cuts!(sepa.cutpool, sepa)
end

function get_cut_from_separating_solution(
    sepa::VPolyhedralSeparator,
    separating_solution::Vector{SCIP.SCIP_Real}
)::Cut
    scip = sepa.scipd
    tableau = sepa.complemented_tableau
    lp_solution = get_solution_vector(tableau)
    lp_solution = get_uncomplemented_vector(lp_solution, tableau)

    separating_solution = undo_projection(sepa.projection, separating_solution)
    separating_solution = get_uncomplemented_vector(separating_solution, tableau)
    b = dot(separating_solution, lp_solution) + 1

    cut_vector, b = convert_standard_inequality_to_general(
        scip, tableau, separating_solution, b
    )

    # We normalize the cut to the form ax <= b
    return Cut(-cut_vector, -b)
end
#=
function solve_separation(
    sepa::VPolyhedralSeparator,
    solution_values::Vector{SCIP.SCIP_Real},
    projection::Projection
)
    scip = sepa.scipd
    lp_solution = get_solution_vector(sepa.complemented_tableau)
    lp_solution = project(projection, lp_solution)
    point_collection = get_points(point_ray_collection)
    point_collection = map(p -> project(projection, p), point_collection)
    ray_collection = get_rays(point_ray_collection)
    ray_collection = map(r -> project(projection, r), ray_collection)
    cut_limit = sepa.cut_limit
    dimension = length(lp_solution)

    separating_lp = Model(HiGHS.Optimizer)
    set_optimizer_attribute(separating_lp, "output_flag", false)

    # Translate the points to nonbasic space
    point_collection = map(p -> p - lp_solution, point_collection)

    # Construct (PRLP0)
    @variable(separating_lp, x[1:dimension] >= 0)
    for point in point_collection
        @constraint(separating_lp, sum(x[i] * point[i] for i in 1:dimension) >= 1)
    end
    for ray in ray_collection
        @constraint(separating_lp, sum(x[i] * ray[i] for i in 1:dimension) >= 0)
    end

    # First check if the LP is feasible by optimizing using all 0s objective
    optimize!(separating_lp)
    if !is_solved_and_feasible(separating_lp)
        error("Separating LP is infeasible")
    end

    # Now optimize with all 1s objective
    @objective(separating_lp, Min, sum(x))
    optimize!(separating_lp)

    if is_solved_and_feasible(separating_lp)
        #push!(cut_collection, value.(x))
    else
        @error ("All ones objective is unbounded")
    end

    # Now Optimize with p as objective
    p_idx = argmin(solution_values)
    p_star = point_collection[p_idx]
    @objective(separating_lp, Min, sum(x[i] * p_star[i] for i in 1:dimension))
    optimize!(separating_lp)

    if is_solved_and_feasible(separating_lp)
        #push!(cut_collection, value.(x))
    else
        @error ("p* as objective is unbounded")
    end

    # Now transition to PRLP=
    a_bar = value.(x)
    @constraint(separating_lp, sum(x[i] * p_star[i] for i in 1:dimension) == 1)

    # points_string = fill("point", length(point_collection))
    # rays_string = fill("ray", length(ray_collection))
    # @info point_score
    # score = [point_score; rays_score]
    # sort = sortperm(score)
    # @info score[sort][1:10]
    # string = [points_string; rays_string]
    # string = string[sort]
    # point_ray = [point_collection; ray_collection]
    # point_ray = point_ray[sort]
    # marker = fill(false, length(point_ray))

    # for i in 1:length(marker)
    #     if marker[i]
    #         continue
    #     end
    #     if string[i] == "point"
    #         if is_EQ(scip, dot(a_bar, point_ray[i]), 1.0)
    #             marker[i] = true
    #         end
    #     else
    #         if is_EQ(scip, dot(a_bar, point_ray[i]), 0.0)
    #             marker[i] = true
    #         end
    #     end
    # end

    # @info string[1:10]
    cnt = 0
    # Ray_collection have been sorted before by orthogonaltiy to the objective
    for r in ray_collection
        @objective(separating_lp, Min, sum(x[i] * r[i] for i in 1:dimension))
        optimize!(separating_lp)

        if cnt < 5
            write_to_file(separating_lp, "separating_lp_$(cnt).lp")
        end

        if cnt % 50 == 0
            @info "Tried $(cnt) rays"
        end
        cnt += 1

        if is_solved_and_feasible(separating_lp)
            sol = value.(x)
            skip = false

            for v in cut_collection
                if is_zero(scip, norm(sol - v))
                    skip = true
                    break
                end
            end

            if skip
                continue
            end

            push!(cut_collection, value.(x))

            # for j in 1:length(marker)
            #     if marker[j]
            #         continue
            #     end
            #     if string[j] == "point"
            #         if is_EQ(scip, dot(a_bar, point_ray[j]), 1.0)
            #             marker[j] = true
            #         end
            #     else
            #         if is_EQ(scip, dot(a_bar, point_ray[j]), 0.0)
            #             marker[j] = true
            #         end
            #     end
            # end

            if length(cut_collection) % 10 == 0
                @info "Number of cut found: $(length(cut_collection))"
            end
        else
            @info solution_summary(separating_lp)
            @warn "Failed to find a cut"
            break
            fail += 1
        end

        if length(cut_collection) >= cut_limit
            break
        end
    end
    @info "Number of cut found: $(length(cut_collection))"
    @info "Disjunctive Lower Bound: $(min(solution_values...))"
    return cut_collection
end
=#