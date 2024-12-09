include("run_single_experiment/Experiment.jl")
include("run_single_experiment/cli_arguments_handlers.jl")
include("run_single_experiment/setup_environment.jl")
include("run_single_experiment/record_result.jl")

experiment = Experiment()
fill_experiment_parameters_from_cli_arguments(experiment)

output_path_relative = get_parameter(experiment, "output_path")
output_path = abspath(output_path_relative)
stdout_path = joinpath(output_path, "stdout.txt")
@info "Output path: $output_path"

setup_scip_parameter(experiment)
load_problem_to_scip(experiment)
#SCIP.@SCIP_CALL SCIP.SCIPpermuteProb(
#    experiment.scip, get_parameter(experiment, "random_seed"), 1, 1, 1, 1, 1
#)
SCIP.@SCIP_CALL SCIP.SCIPsolve(experiment.scip)
write_output(experiment)
