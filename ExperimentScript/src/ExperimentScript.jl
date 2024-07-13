module ExperimentScript
include("UtilityStructures.jl")
include("commandline_handler.jl")
include("separator_wrappers.jl")
include("scip_settings.jl")
include("setup_environment.jl")
include("load_problem_data.jl")
include("postprocessing.jl")

function run_experiment()
    execution_parameters = get_execution_parameters()
    experiment_store = setup_environment(execution_parameters)
    load_problem_data(experiment_store)
    SCIP.@SCIP_CALL SCIP.SCIPsolve(experiment_store.scip)
    do_postprocessing(experiment_store)
end
end # module ExperimentScript
