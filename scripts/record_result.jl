function record_result(experiment_store::ExperimentStore)
    write_output(experiment_store)
end

function write_output(experiment_store::ExperimentStore)
    output = open(joinpath(experiment_store.result_path, "results.txt"), "w")

    println(output, "ExperimentMode: $(experiment_store.separator)")
    println(output, "Instance: $(experiment_store.instance)")
    println(output, "OutputDirectory: $(experiment_store.result_path)")
    println(output, "SCIPStatus: $(SCIP.SCIPgetStatus(experiment_store.scip))")

    println(
        output,
        "ReferenceObjective: $(round(experiment_store.reference_objective,sigdigits = 6))"
    )

    initialdualbound = SCIP.SCIPgetFirstLPDualboundRoot(experiment_store.scip)
    println(output, "FirstLPObj: $(round(initialdualbound,sigdigits = 6))")

    initialgap = abs(initialdualbound - experiment_store.reference_objective)
    println(output, "InitialGap: $(round(initialgap,sigdigits = 6))")

    finaldualbound = SCIP.SCIPgetDualboundRoot(experiment_store.scip)
    println(output, "FinalLPObj: $(round(finaldualbound,sigdigits = 6))")

    finalgap = abs(finaldualbound - experiment_store.reference_objective)
    println(output, "FinalGap: $(round(finalgap, sigdigits=6 ))")

    gapclosed = 1 - finalgap / initialgap
    println(output, "GapClosed: $(round(gapclosed, sigdigits=6))")

    check = false
    if (SCIP.SCIPgetObjsense(experiment_store.scip) == SCIP.SCIP_OBJSENSE_MINIMIZE)
        check = experiment_store.reference_objective > finaldualbound
    elseif (SCIP.SCIPgetObjsense(experiment_store.scip) == SCIP.SCIP_OBJSENSE_MAXIMIZE)
        check = experiment_store.reference_objective < finaldualbound
    end
    println(output, "Check: $(check)")

    @info "GapClosed: $(round(gapclosed, sigdigits=3)*100)%"
    if check
        @info "Check: OK"
    else
        @error "Check: Failed"
    end
    write_solving_statistic(experiment_store)
end

function write_solving_statistic(experiment_store::ExperimentStore)
    solv_stats = open(joinpath(experiment_store.result_path, "solving_statistic.txt"), "w")
    solv_stats_ptr = Libc.FILE(solv_stats)
    SCIP.SCIPprintStatistics(experiment_store.scip, solv_stats_ptr)
end