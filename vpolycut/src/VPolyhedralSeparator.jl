using SCIP
using JuMP
using LinearAlgebra
import MathOptInterface as MOI
using HiGHS

"""
VPolyhedral Cut Separator

Implementation of Algorithm 1 from 
Balas, Egon, and Aleksandr M. Kazachkov. "V-polyhedral disjunctive cuts." arXiv preprint arXiv:2207.13619 (2022).
Disjunction are obtained from partial branch and bound trees

# Required Parameters
- scipd::SCIP.SCIPData Reference to the SCIPData object

# Optional Parameters
- n_leaves::Int Number of leaves in the disjunction
- cut_limit::Int Number of cuts to generate (-1 if no limit, -2 if limit is the number of fractional variable)
"""
@kwdef mutable struct VPolyhedralSeparator <: SCIP.AbstractSeparator
    "Pointer to SCIP"
    scipd::SCIP.SCIPData
    "Number of leaves in the disjunction"
    n_leaves::Int = 2
    "Number of times the separator is called"
    called::Int = 0
    "Cut Limit"
    cut_limit::Int = -2
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
    "LP Solution"
    lp_sol::Union{Nothing,Vector{SCIP.SCIP_Real}} = nothing
end

# Include Helper
function include_vpolyhedral_sepa(scip::SCIP.SCIPData; n_leaves=2, cut_limit=-2)
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
    # If cut limit is -1 or -2 convert them to the actual limit 
    println(sepa.cut_limit)
    sepa.cut_limit = get_cut_limit(sepa)
    println("Cut limit is $(sepa.cut_limit)")
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
end

function get_cut_limit(sepa::VPolyhedralSeparator)
    if sepa.cut_limit == -1
        return typemax(Int)
    elseif sepa.cut_limit == -2
        return SCIP.SCIPgetNLPBranchCands(sepa.scipd)
    else
        return sepa.cut_limit
    end
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
        println("Generated $(length(sepa.cutpool)) cuts. Cutlimit is $(sepa.cut_limit)")
        if length(sepa.cutpool) >= sepa.cut_limit
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