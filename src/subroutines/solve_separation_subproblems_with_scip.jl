using SCIP

function solve_separation_subproblems(
    scip::SCIP.SCIPData, points::Vector{Point}, rays::Vector{Ray}, p_star::Point,
    statistic::Vector{Any}, cut_limit::Int64, time_limit::Float64
)
    # Setup aliases
    start_time = time()

    # Get current LP Solution projected onto the non-basic variables.
    # This will be the origin point 
    problem_dimension = length(p_star)
    solution_pool = []

    points_zeroed = map(points) do point
        return [!is_zero(scip, p) ? p : 0.0 for p in point]
    end

    rays_zeroed = map(rays) do ray
        return [!is_zero(scip, p) ? p : 0.0 for p in get_coefficients(ray)]
    end

    # Construct the PRLP
    @debug "Constructing PRLP"
    prlp = PRLP(problem_dimension)

    for point in points_zeroed
        PRLPaddPoint(prlp, point)
    end
    for ray in rays_zeroed
        PRLPaddRay(prlp, ray)
    end

    PRLPconstructLP(prlp)
    @debug "PRLP Constructed"

    @debug "Checking PRLP Feasibility"
    PRLPsetTimeLimit(prlp, 300.0)

    feasibility_check = false
    # We iterate over the algorithms to check which one find a feasible solution
    for algorithm in [PRIMAL_SIMPLEX, DUAL_SIMPLEX, BARRIER]
        PRLPsetSolvingAlgorithm(prlp, algorithm)
        @debug "Trying Algorithm $(algorithm)"

        feasibility_check = try_objective(
            solution_pool, statistic, prlp, zeros(SCIP.SCIP_Real, problem_dimension),
            "feasibility"
        )
        if feasibility_check
            @debug "Success with $(algorithm)"
            break
        end
    end
    if !feasibility_check
        PRLPdestroy(prlp)
        throw(FailedToProvePRLPFeasibility())
    end

    # All Ones
    PRLPsetTimeLimit(prlp, 30.0)
    try_objective(
        solution_pool, statistic, prlp, ones(SCIP.SCIP_Real, problem_dimension), "all_ones"
    )
    # We can live if all ones is not optimized to optimality 

    # Pstar
    p_run = try_objective(solution_pool, statistic, prlp, p_star, "p_star")
    # Unless it is a barrier solve, only proceed if objective is 1.
    if !p_run && PRLPgetSolvingAlgorithm(prlp) != BARRIER &&
        is_EQ(scip, PRLPgetObjValue(prlp), 1.0)
        PRLPdestroy(prlp)
        throw(PStarNotTight())
    end

    # Check if tightened Pstar is feasible
    @debug "Trying Pstar Tightening"
    PRLPsetTimeLimit(prlp, 300.0)
    p_star_zeroed = [!is_zero(scip, p) ? p : 0.0 for p in p_star]
    PRLPtighten(prlp, p_star_zeroed)
    tight_try = try_objective(
        solution_pool, statistic, prlp, zeros(SCIP.SCIP_Real, problem_dimension),
        "prlp=_feasibility"
    )
    if !tight_try
        PRLPdestroy(prlp)
        throw(PStarInfeasible())
    end

    @debug "Pstar Tightening Successful"
    PRLPsetTimeLimit(prlp, 30.0)
    a_bar = PRLPgetSolution(prlp)
    r_bar = filter(ray -> !is_zero(scip, dot(a_bar, ray)), rays)
    sort!(r_bar; by=ray -> abs(get_obj(get_generating_variable(ray))))
    objective_tried = 0

    for ray in r_bar
        r_run = try_objective(
            solution_pool, statistic, prlp, get_coefficients(ray),
            "prlp=_$(objective_tried)"
        )
        objective_tried += 1
        if r_run
            @debug "Generated $(length(solution_pool)) cuts. Cutlimit is $(cut_limit)"
        end
        if length(solution_pool) >= cut_limit
            @goto cleanup
        end
        elapsed_time = time() - start_time
        if elapsed_time > 900
            @goto cleanup
        end
    end

    @label cleanup
    PRLPdestroy(prlp)
    return solution_pool
end

function get_cut_from_separating_solution(
    scip::SCIP.SCIPData, tableau::Tableau, nonbasicspace::NonBasicSpace,
    separating_solution::Vector{SCIP.SCIP_Real}
)::Cut
    lp_solution = get_solution_vector(tableau)

    separating_solution = revert_point_to_original_space(nonbasicspace, separating_solution)
    b = dot(separating_solution, lp_solution) + 1

    cut_vector, b = convert_standard_inequality_to_general(
        scip, tableau, separating_solution, b
    )

    # We normalize the cut to the form ax <= b
    return Cut(-cut_vector, -b)
end

function try_objective(
    solution_pool::Vector{Any}, statistic::Vector{Any}, prlp::PRLP,
    objective::Vector{SCIP.SCIP_Real}, label=""
)
    PRLPsetObjective(prlp, objective)
    PRLPsolve(prlp)
    push!(statistic, get_solve_stat(prlp, label))
    if PRLPisSolutionAvailable(prlp)
        push!(solution_pool, PRLPgetSolution(prlp))
        return true
    end
    return false
end

function get_solve_stat(prlp::PRLP, obj_name::String)
    # Since SCIP LPI doesn't support solve time, this needed to be passed manually
    return Dict(
        "solve_algorithm" => PRLPgetSolvingAlgorithm(prlp),
        "solve_time" => PRLPgetLastSolveTime(prlp),
        "simplex_iterations" => PRLPgetLastSimplexIterations(prlp),
        "primal_status" => PRLPisSolutionAvailable(prlp),
        "termination_status" => PRLPgetLastTerminationStatus(prlp),
        "objective_name" => obj_name
    )
end
