include("cli_arguments.jl")
include("separator_wrappers.jl")
include("scip_settings.jl")
include("setup_environment.jl")
include("load_problem_data.jl")
include("record_result.jl")

function run_experiment()
    execution_parameters = get_execution_parameters()
    experiment_store = setup_environment(execution_parameters)
    load_problem_data(experiment_store)
    SCIP.@SCIP_CALL SCIP.SCIPsolve(experiment_store.scip)
    record_result(experiment_store)
end

run_experiment()