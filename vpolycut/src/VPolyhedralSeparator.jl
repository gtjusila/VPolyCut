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
"""
@kwdef mutable struct VPolyhedralSeparator <: SCIP.AbstractSeparator
    "Pointer to SCIP"
    scipd::SCIP.SCIPData
    "Number of times the separator is called"
    called::Int = 0
    "Number of leaves in the disjunction"
    n_leaves::Int = 2
    "Is LP solution separated?"
    separated::Bool = false
end

"""
The the point ray collection and solution value from a disjunctive term
The point ray collection will be given in the complemented space therefore the 
original root_tableau is required.

# Parameters
- scip::SCIP.SCIPData Reference to the SCIPData object
- root_tableau::ComplementedTableau The root tableau of the problem
- path::Vector{Node} The path from the root to the disjunctive term

# Returns
- points_collection::Vector{Point} A corner polyhedron vertex 
- rays_collection::Vector{Ray} The collection of corner polyhedron rays
- solution_values::Vector{SCIP.SCIP_Real} The solution value at the vertex
"""
function get_disjunctive_term_information(
    scip::SCIP.SCIPData,
    root_tableau::ComplementedTableau,
    path::Vector{Node}
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
    #@info "Tableau" tableau.tableau_matrix [
    #    get_basis_status(get_var_from_column(tableau, i)) for i in 1:4
    #]

    # Complement the same columns as the root tableau
    complemented_columns = get_complemented_columns(root_tableau)
    for i in complemented_columns
        var = get_var_from_column(tableau, i)
        complement_column(tableau, var)
    end
    #@info "Complemented Tableau" tableau.tableau_matrix
    corner = construct_corner_polyhedron(tableau)
    solution = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    @assert get_nvars(tableau) == get_nvars(root_tableau)

    # Leave Probing mode
    SCIP.SCIPendProbing(scip)
    return [get_lp_sol(corner)], get_lp_rays(corner), [solution]
end

function get_cuts(
    scip::SCIP.SCIPData,
    lp_solution::Point,
    point_collection::Vector{Point},
    ray_collection::Vector{Ray},
    solution_values::Vector{SCIP.SCIP_Real},
    cut_limit::Int
)
    dimension = length(lp_solution)
    @info "Dimension: $dimension"
    separating_lp = Model(HiGHS.Optimizer)
    set_optimizer_attribute(separating_lp, "output_flag", false)
    cut_collection = Point[]

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

    # Now Optimize with p* as objective
    min_ind = argmin(solution_values)
    p_star = point_collection[min_ind]
    @objective(separating_lp, Min, sum(x[i] * p_star[i] for i in 1:dimension))
    optimize!(separating_lp)

    if is_solved_and_feasible(separating_lp)
        #push!(cut_collection, value.(x))
    else
        @error ("p* as objective is unbounded")
    end

    # Now transition to PRLP=
    a_bar = value.(x)
    @constraint(separating_lp, sum(x[i] * p_star[i] for i in 1:dimension) <= 1)

    # Ray_collection have been sorted before by orthogonaltiy to the objective
    cnt = 0
    r_bar = filter(r -> is_non_zero(scip, dot(a_bar, r)), ray_collection)
    r = popfirst!(r_bar)
    while true
        @objective(separating_lp, Min, sum(x[i] * r[i] for i in 1:dimension))
        optimize!(separating_lp)
        cnt += 1
        if cnt < 5
            write_to_file(separating_lp, "separating_lp_$(cnt).lp")
        end

        if is_solved_and_feasible(separating_lp)
            push!(cut_collection, value.(x))
            sol = value.(x)
            @info sol[1:5] (dot(sol, a_bar) / (norm(sol) * norm(a_bar))) (norm(sol - a_bar))
            a_bar = sol
            @info "Cut Found objective: $(objective_value(separating_lp)) Pdot = $(sum(sol[i]*p_star[i] for i in 1:dimension)) r_bar length=$(length(r_bar))"
            if length(cut_collection) % 10 == 0
                @info "Number of cut found: $(length(cut_collection))"
            end
        else
            @info solution_summary(separating_lp)
            @warn "Failed to find a cut"
            fail += 1
        end
        r_bar = filter(p -> (abs(dot(r, p)) < 0.1), r_bar)
        if isempty(r_bar)
            break
        end
        r = popfirst!(r_bar)
        if length(cut_collection) >= cut_limit
            break
        end
    end

    @info "Number of cut found: $(length(cut_collection))"
    @info "Disjunctive Lower Bound: $(solution_values[min_ind])"
    return cut_collection
end

function vpolyhedralcut_separation(sepa::VPolyhedralSeparator)
    #
    # Preparation
    #
    scip = sepa.scipd
    # For experimental purpose limit the number of cuts to the number of branching candidate
    cut_limit = Int(SCIP.SCIPgetNLPBranchCands(scip))

    # Step 0: Get complemented tableau
    original_tableau = construct_tableau_with_constraint_matrix(scip)
    complemented_tableau = ComplementedTableau(original_tableau)
    lp_sol = get_solution_vector(complemented_tableau)
    objective_direction = get_objective_direction(complemented_tableau)

    # Step 1: Get Disjunction
    @info "Starting Branch And Bound to generate disjunction"
    branchandbound = BranchAndBound(
        scip; max_leaves=64
    )
    execute_branchandbound(branchandbound)
    disjunctions = get_leaves(branchandbound)
    @info "Branch and bound completed"

    #
    # Main Algorithm 
    #

    # Step 1: Get the points and rays from the disjunctions
    @info "Collecting Point Ray Collection"
    points_collection = Point[]
    rays_collection = Ray[]
    solution_values = SCIP.SCIP_Real[]
    for term in disjunctions
        path = get_path(term) # Get all actions from root to leaf
        point, rays, solution_value = get_disjunctive_term_information(
            scip, complemented_tableau, path
        )
        rays = filter(r -> is_non_zero(scip, norm(r)), rays)
        append!(points_collection, point)
        append!(rays_collection, rays)
        append!(solution_values, solution_value)
    end

    # Sort Ray By Orthogonality to the objective if objective is not the zero vector
    if is_non_zero(scip, norm(objective_direction))
        sort!(
            rays_collection;
            by=r -> abs(dot(objective_direction, r)),
            rev=false)
    end
    @info "Point Ray Collection Collected"

    # Step 2: Project the lp_solution, points and rays to the non-basic space
    projection = create_projection_to_nonbasic_space(complemented_tableau)
    projected_lp_sol = project_point(projection, lp_sol)
    projected_points = [project_point(projection, point) for point in points_collection]
    projected_rays = [project_point(projection, ray) for ray in rays_collection]

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

function include_vpolyhedral_sepa(scip::SCIP.SCIPData; n_leaves=2)
    sepa = VPolyhedralSeparator(; scipd=scip, n_leaves=n_leaves)
    SCIP.include_sepa(
        scip.scip[], scip.sepas, sepa; priority=9999, freq=0, usessubscip=true
    )
end
