using SCIP
using JSON
using VPolyhedralCut

function write_output(experiment::Experiment)
    output_path_relative = get_parameter(experiment, "output_path")
    output_path = abspath(output_path_relative)

    # Experiment Result
    result = Dict{String,Any}()
    result["separator"] = experiment.parameters.separator
    result["instance"] = experiment.parameters.instance_path
    result["scip_status"] = SCIP.SCIPgetStatus(experiment.scip)
    result["initial_lp_obj"] = SCIP.SCIPgetFirstLPDualboundRoot(experiment.scip)
    result["final_lp_obj"] = SCIP.SCIPgetDualboundRoot(experiment.scip)

    result_path = joinpath(output_path, "results.json")
    open(result_path, "w") do io
        JSON.print(io, result, 4)
    end

    # SCIP Solving Statistic
    solving_statistic_path = joinpath(output_path, "scip_solving_statistic.txt")
    open(solving_statistic_path, "w") do io
        solv_stats_ptr = Libc.FILE(io)
        SCIP.SCIPprintStatistics(experiment.scip, solv_stats_ptr)
    end
end

#=
function write_solving_statistic(experiment_store::ExperimentStore)
    solv_stats_path = open(
        joinpath(experiment_store.result_path, "solving_statistic.txt"), "w"
    )
    solv_stats_ptr = Libc.FILE(solv_stats_path)
    SCIP.SCIPprintStatistics(experiment_store.scip, solv_stats_ptr)
end
=#