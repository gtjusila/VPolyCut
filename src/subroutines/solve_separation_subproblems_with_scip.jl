using SCIP

function solve_separation_subproblems(sepa::VPCSeparator)
    # Setup aliases
    scip = sepa.scipd

    # Get current LP Solution projected onto the non-basic variables.
    # This will be the origin point 
    current_lp_solution = get_solution_vector(sepa.complemented_tableau)
    projected_current_lp_solution = project(sepa.projection, current_lp_solution)
    problem_dimension = length(projected_current_lp_solution)

    # First, translate the points so that the lp_solution is at the origin 
    projected_point_collection = get_points(sepa.point_ray_collection)
    point_collection_in_nonbasic_space = map(projected_point_collection) do corner_point
        return CornerPoint(
            get_point(corner_point) - projected_current_lp_solution,
            get_objective_value(corner_point),
            get_orig_objective_value(corner_point)
        )
    end
    ray_collection = get_rays(sepa.point_ray_collection)

    # Construct the PRLP
    @debug "Constructing PRLP"
    prlp = PRLP(problem_dimension)

    for point in point_collection_in_nonbasic_space
        PRLPaddPoint(prlp, get_point(point))
    end
    for ray in ray_collection
        PRLPaddRay(prlp, get_coefficients(ray))
    end
    PRLPconstructLP(prlp)
    @debug "PRLP Constructed"

    @debug "Checking PRLP Feasibility"
    PRLPsetTimeLimit(prlp, 300.0)
    feasibility_check = false
    for algorithm in [PRIMAL_SIMPLEX, BARRIER]
        PRLPsetSolvingAlgorithm(prlp, algorithm)
        @debug "Trying Algorithm $(algorithm)"
        feasibility_check = try_objective(
            sepa, prlp, zeros(SCIP.SCIP_Real, problem_dimension), "feasibility"
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
        sepa, prlp, ones(SCIP.SCIP_Real, problem_dimension), "all_ones"
    )
    # We can live if all ones is not optimized to optimality 

    # Pstar
    p_star = argmin(point -> get_objective_value(point), point_collection_in_nonbasic_space)
    p_star = get_point(p_star)
    p_run = try_objective(sepa, prlp, p_star, "p_star")
    # Unless it is a barrier solve, only proceed if objective is 1.
    if !p_run && PRLPgetSolvingAlgorithm(prlp) != BARRIER &&
        is_EQ(scip, PRLPgetObjValue(prlp), 1.0)
        PRLPdestroy(prlp)
        throw(PStarNotTight())
    end

    # Check if tightened Pstar is feasible
    @debug "Trying Pstar Tightening"
    PRLPsetTimeLimit(prlp, 300.0)
    PRLPtighten(prlp, p_star)
    tight_try = try_objective(
        sepa, prlp, zeros(SCIP.SCIP_Real, problem_dimension), "prlp=_feasibility"
    )
    if !tight_try
        PRLPdestroy(prlp)
        throw(PStarInfeasible())
    end

    @debug "Pstar Tightening Successful"
    PRLPsetTimeLimit(prlp, 30.0)
    a_bar = PRLPgetSolution(prlp)
    r_bar = filter(ray -> !is_zero(scip, dot(a_bar, ray)), ray_collection)
    sort!(r_bar; by=ray -> abs(get_obj(get_generating_variable(ray))))
    objective_tried = 0

    for ray in r_bar
        r_run = try_objective(sepa, prlp, get_coefficients(ray), "prlp=_$(objective_tried)")
        objective_tried += 1
        if r_run
            @debug "Generated $(length(sepa.cutpool)) cuts. Cutlimit is $(sepa.parameters.cut_limit)"
        end
        if length(sepa.cutpool) >= sepa.parameters.cut_limit
            @goto cleanup
        end
        elapsed_time = time() - sepa.start_time
        if elapsed_time > sepa.parameters.time_limit
            @goto cleanup
        end
    end

    @label cleanup
    PRLPdestroy(prlp)
    add_all_cuts!(sepa.cutpool, sepa)
    if length(sepa.cutpool) >= 0
        sepa.separated = true
    end
end

function get_cut_from_separating_solution(
    sepa::VPCSeparator,
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

function try_objective(
    sepa::VPCSeparator, prlp::PRLP, objective::Vector{SCIP.SCIP_Real}, label=""
)
    PRLPsetObjective(prlp, objective)
    PRLPsolve(prlp)
    push!(sepa.prlp_solves, get_solve_stat(prlp, label))
    if PRLPisSolutionAvailable(prlp)
        push!(sepa.cutpool, get_cut_from_separating_solution(sepa, PRLPgetSolution(prlp)))
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
