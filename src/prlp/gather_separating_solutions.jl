function gather_separating_solutions(
    sepa::VPCSeparator
)
    shared_data = sepa.shared_data
    statistic = sepa.statistics

    # Data that is needed from VPCSeparator
    prlp = sepa.shared_data.prlp
    point_ray_collection = sepa.shared_data.point_ray_collection
    cut_limit = sepa.parameters.cut_limit
    time_limit = sepa.parameters.time_limit
    start_time = sepa.shared_data.start_time

    @debug "Remaining Time for Separation: $(time_limit - (time() - start_time))"

    separating_solutions = []
    objective_pool = ObjectivePool(prlp, point_ray_collection)
    total_time = 0.0
    consecutive_fail = 0
    disjunctive_gap = shared_data.disjunctive_lower_bound - shared_data.lp_obj
    @info "Disjunctive Gap $(disjunctive_gap)"

    SCIP.@SCIP_CALL SCIP.SCIPstartDive(shared_data.scipd)
    for (counter, objective) in enumerate(objective_pool)
        @debug "Trying objective $(objective.label)"
        objective_solution_timed = @timed PRLPtryObjective(prlp, objective)
        objective_solution = objective_solution_timed.value
        total_time += objective_solution_timed.time
        if !isnothing(objective_solution)
            consecutive_fail = 0
            push!(separating_solutions, objective_solution)
            cut = get_cut_from_separating_solution(
                objective_solution, sepa.shared_data.nonbasic_space
            )
            add_diving_row!(
                shared_data.scipd,
                sepa,
                get_coefficients(cut),
                shared_data.nonbasic_space.variable_pointers,
                get_rhs(cut)
            )

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
            percent_disjunctive_gap_closed =
                ((lp_obj - shared_data.lp_obj) / disjunctive_gap) * 100

            push!(
                statistic.prlp_percent_disjunctive_gap_closed_history,
                percent_disjunctive_gap_closed
            )
            @debug "Found cut. Total cut found so far: $(length(separating_solutions)). Cut limit is $(cut_limit). percent disjunctive gap closed $(percent_disjunctive_gap_closed)"
        else
            consecutive_fail += 1
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
        # Objective limit reached
        if counter >= 2 * cut_limit
            break
        end
        if consecutive_fail >= 10
            # Mark but still fail softly
            sepa.termination_status = CONSECUTIVE_FAIL_LIMIT_REACHED
            break
        end
    end
    SCIP.@SCIP_CALL SCIP.SCIPendDive(shared_data.scipd)
    PRLPfreeBasis(prlp)
    return separating_solutions
end