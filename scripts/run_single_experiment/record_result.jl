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
    if result["separator"] == "vpc"
        result["sepa_termination_message"] = experiment.vpcsepa.termination_message
        result["number_of_cuts"] = length(experiment.vpcsepa.cutpool)
        result["disjunctive_lower_bound"] = experiment.vpcsepa.disjunctive_lower_bound
        result["n_fractional_variables"] = experiment.vpcsepa.n_fractional_variables
        result["prlp_solves"] = experiment.vpcsepa.prlp_solves
        result["cbar_test"] = experiment.vpcsepa.cbar_test
    else
        result["sepa_termination_message"] = ""
        result["number_of_cuts"] = -1
        result["disjunctive_lower_bound"] = -1.0
        result["n_fractional_variables"] = -1
        result["prlp_solves"] = []
        result["cbar_test"] = false
    end
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