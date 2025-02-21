function gather_separating_solutions(
    prlp::PRLP,
    point_ray_collection::PointRayCollection;
    cut_limit::Int = typemax(Int),
    time_limit::Float64 = typemax(Float64),
    start_time::Float64 = time()
)
    @debug "Remaining Time for Separation: $(time_limit - (time() - start_time))"
    separating_solutions = []
    objective_pool = ObjectivePool(prlp, point_ray_collection)
    total_time = 0.0
    consecutive_fail = 0
    for (counter, objective) in enumerate(objective_pool)
        @debug "Trying objective $(objective.label)"
        objective_solution_timed = @timed PRLPtryObjective(prlp, objective)
        objective_solution = objective_solution_timed.value
        total_time += objective_solution_timed.time
        if !isnothing(objective_solution)
            consecutive_fail = 0
            push!(separating_solutions, objective_solution)
            @debug "Found cut. Total cut found so far: $(length(separating_solutions)). Cut limit is $(cut_limit)"
        else
            consecutive_fail += 1
        end
        # Time limit exceeded
        if time_limit < time() - start_time
            sepa.termination_status = TIME_LIMIT_EXCEEDED
            break
        end
        # Cut limit reached
        if length(separating_solutions) >= cut_limit
            break
        end
        # Objective limit reached
        if counter >= 2 * cut_limit
            break
        end
        if consecutive_fail >= 10
            sepa.termination_status = CONSECUTIVE_FAIL_LIMIT_REACHED
            break
        end
    end
    PRLPfreeBasis(prlp)
    return separating_solutions
end