function gather_separating_solutions(
    sepa::VPCSeparator
)
    shared_data = sepa.shared_data
    statistic = sepa.statistics

    # Data that is needed from VPCSeparator
    prlp = sepa.shared_data.prlp
    point_ray_collection = sepa.shared_data.point_ray_collection
    cut_limit = sepa.parameters.cut_limit
    max_consecutive_fail = sepa.parameters.prlp_max_consecutive_fail
    min_increase = sepa.parameters.prlp_min_increase_non_stagnating
    max_stagnating_rounds = sepa.parameters.prlp_max_consecutive_stagnation
    time_limit = sepa.parameters.prlp_time_limit
    start_time = time()

    @debug "Remaining Time for Separation: $(time_limit - (time() - start_time))"

    # Initialize Variables
    separating_solutions = []
    objective_pool = ObjectivePool(
        prlp, point_ray_collection, Vector{Point}([])
    )
    consecutive_fail = 0
    disjunctive_gap = shared_data.disjunctive_lower_bound - shared_data.lp_obj
    gap_closed_history::Vector{SCIP.SCIP_Real} = []
    @info "Disjunctive Gap $(disjunctive_gap)"

    # Do everything in dive mode
    SCIP.@SCIP_CALL SCIP.SCIPstartDive(shared_data.scipd)

    # Loop to gather separating solutions   
    for (counter, objective) in enumerate(objective_pool)
        @debug "Trying objective $(objective.label)"
        objective_solution_timed = @timed PRLPtryObjective(prlp, objective)
        objective_solution = objective_solution_timed.value

        if !isnothing(objective_solution)
            consecutive_fail = 0
            push!(separating_solutions, objective_solution)
            cut = get_cut_from_separating_solution(
                objective_solution, sepa.shared_data.nonbasic_space, prlp.beta
            )
            add_diving_row!(
                shared_data.scipd,
                sepa,
                get_coefficients(cut),
                shared_data.nonbasic_space.variable_pointers,
                get_rhs(cut)
            )

            @debug "Found cut. Total cut found so far: $(length(separating_solutions)). Cut limit is $(cut_limit)."
        else
            consecutive_fail += 1
        end

        # Get new LP relaxation objective
        lp_obj = get_new_lp_relaxation_objective(sepa)
        relative_disjunctive_gap_closed = (lp_obj - shared_data.lp_obj) / disjunctive_gap
        push!(
            gap_closed_history,
            relative_disjunctive_gap_closed
        )
        @info "Relative Disjunctive Gap Closed: $(relative_disjunctive_gap_closed)"

        # Check if we are stagnating
        if (counter > max_stagnating_rounds)
            closed =
                relative_disjunctive_gap_closed -
                gap_closed_history[counter - max_stagnating_rounds]
            if is_LT(closed, min_increase)
                @info "Stagnating. Stopping separation."
                break
            end
        end
        # Time limit exceeded
        if time_limit < time() - start_time
            # Mark as time limit reached but we still fail softly to allow cuts to be added 
            sepa.termination_status = TIME_LIMIT_EXCEEDED_PRLP
            break
        end
        # Cut limit reached
        if length(separating_solutions) >= cut_limit
            break
        end
        if consecutive_fail >= max_consecutive_fail
            # Mark but still fail softly
            sepa.termination_status = CONSECUTIVE_FAIL_LIMIT_REACHED
            break
        end
    end
    statistic.prlp_percent_disjunctive_gap_closed_history = gap_closed_history * 100
    SCIP.@SCIP_CALL SCIP.SCIPendDive(shared_data.scipd)
    PRLPfreeBasis(prlp)
    return separating_solutions
end

function get_new_lp_relaxation_objective(sepa::VPCSeparator)
    shared_data = sepa.shared_data
    statistic = sepa.statistics

    lp_error = Ref{SCIP.SCIP_Bool}(false)
    cutoff = Ref{SCIP.SCIP_Bool}(false)
    SCIP.SCIPsolveDiveLP(shared_data.scipd, -1, lp_error, cutoff)

    if is_true(lp_error[])
        throw(LPError())
    end

    if is_true(cutoff[])
        throw(LPError())
    end

    lp_obj = SCIP.SCIPgetSolOrigObj(shared_data.scipd, C_NULL)
    return lp_obj
end