using ArgParse
using UUIDs
function fill_experiment_parameters_from_cli_arguments(experiment::Experiment)
    cli_arguments = setup_cli_arguments()

    if isnothing(cli_arguments["output"])
        cli_arguments["output"] = setup_temporary_output_directory()
    end

    set_parameter(experiment, "output_path", cli_arguments["output"])
    set_parameter(experiment, "separator", cli_arguments["mode"])
    set_parameter(experiment, "instance_path", cli_arguments["instance_path"])
    set_parameter(experiment, "presolving", !cli_arguments["easy"])
    set_parameter(experiment, "propagation", !cli_arguments["easy"])
    set_parameter(experiment, "time_limit", parse(Int64, cli_arguments["time_limit"]))
    set_parameter(experiment, "zeroing_heuristic", cli_arguments["zeroing_heuristic"])
end

function setup_cli_arguments()
    settings = ArgParseSettings()
    @add_arg_table settings begin
        "--instance_path", "-i"
        help = "Path to instance"
        required = true
        "--mode", "-m"
        help = "The mode of experiment to be run. Modes available: gomory, vpoly"
        required = true
        "--easy", "-e"
        help = "Disable Presolving and Propagation for easy instances"
        action = :store_true
        "--output", "-o"
        help = "Output directory for the results"
        "--time_limit", "-t"
        help = "Time limit for the experiment"
        default = "7200"
        "--zeroing_heuristic", "-z"
        help = "Enable Zeroing Heuristic"
        action = :store_true
    end
    return parse_args(settings)
end

function setup_temporary_output_directory()
    tmp = joinpath(pwd(), "temp")
    if !isdir(tmp)
        mkdir(tmp)
    end

    random_string = string(uuid4())
    tmp_path = joinpath(tmp, random_string)
    if !isdir(tmp_path)
        mkdir(tmp_path)
    end

    return tmp_path
end