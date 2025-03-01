using DataFrames
using CSV
using Mustache
using TOML
include("utils.jl")

println("Setup experiment")

instance_list = prompt_user(;
    message = "Instance List",
    validation = (x) -> isfile(abspath(x)),
    error_message = "Invalid Path.",
    default = "experiment_data/miplibbench/instances_list.txt"
)

instance_dir = prompt_user(;
    message = "Instance List",
    validation = (x) -> isdir(abspath(x)),
    error_message = "Invalid Path.",
    default = "experiment_data/miplibbench"
)

mode = prompt_user(;
    message = "Mode",
    validation = (x) -> x in ["vpc", "gomory"],
    error_message = "Invalid Mode.",
    default = "vpc"
)

# If experiment_runs is not there create it first.
runs_path = abspath(joinpath(dirname(@__FILE__), "..", "experiment_runs"))
if !isdir(runs_path)
    @info "Creating experiment_runs folder"
    mkdir(runs_path)
end

experiment_path = prompt_user(;
    message = "Experiment Name:",
    validation = (x) ->
        ((x != "") && !isdir(joinpath(runs_path, x))),
    error_message = "Name must not be empty.",
    default = ""
)
experiment_path = joinpath(runs_path, experiment_path)
mkdir(experiment_path)

### Done with basic parameter, start reading input ###
instances = readlines(instance_list)

if mode == "gomory"
    ### Creating Config TSV ###
    instances_path = [joinpath(instance_dir, instance * ".mps") for instance in instances]
    output_path = [joinpath(experiment_path, instance) for instance in instances]
    for path in output_path
        mkdir(path)
    end
    config_dataframe = DataFrame(;
        id = 1:length(instances),
        instances_path = instances_path,
        output_path = output_path
    )
    run_settings_file = joinpath(experiment_path, "experiment_list.tsv")
    CSV.write(run_settings_file, config_dataframe; delim = '\t')

    for path in instances_path
        if (!isfile(path))
            @warn "Possibly missing instance $(path)"
        end
    end

    # Read the bash template file and write to the experiment folder
    template_file = joinpath(@__DIR__, "gomory_template.sh")
    data = Dict(
        "N" => length(instances),
        "EXPERIMENT_PATH" => experiment_path
    )
    bash_template = read(template_file, String)
    bash_script = Mustache.render(bash_template, data)
    output_file = joinpath(experiment_path, "job_script.sh")
    open(output_file, "w") do f
        write(f, bash_script)
    end
    run(`sbatch $output_file`)

elseif mode == "vpc"
    ### Creating Config TSV ###
    instances_path = [joinpath(instance_dir, instance * ".mps") for instance in instances]
    output_path = [joinpath(experiment_path, instance) for instance in instances]
    for path in output_path
        mkdir(path)
    end
    config_dataframe = DataFrame(;
        id = 1:length(instances),
        instances_path = instances_path,
        output_path = output_path
    )
    run_settings_file = joinpath(experiment_path, "experiment_list.tsv")
    CSV.write(run_settings_file, config_dataframe; delim = '\t')

    for path in instances_path
        if (!isfile(path))
            @warn "Possibly missing instance $(path)"
        end
    end

    # Create config file
    vpc_config = Dict()
    vpc_config["n_leaves"] = prompt_user(;
        message = "Number of Leaves",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        default = "64"
    )
    vpc_config["n_leaves"] = parse(Int, vpc_config["n_leaves"])
    vpc_config["prlp_solve_method"] = prompt_user(;
        message = "PRLP Solve Method (1: PRIMAL SIMPLEX, 2:DUAL SIMPLEX, 3: BARRIER)",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        default = "1"
    )
    vpc_config["prlp_solve_method"] = parse(Int, vpc_config["prlp_solve_method"])
    vpc_config["prlp_allow_warm_start"] = prompt_user(;
        message = "Allow PRLP warm start (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        default = "true"
    )
    vpc_config["prlp_allow_warm_start"] = (vpc_config["prlp_allow_warm_start"] == "true")
    vpc_config["disable_scip_cuts"] = prompt_user(;
        message = "Disable SCIP Cuts (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        default = "true"
    )
    vpc_config["disable_scip_cuts"] = (vpc_config["disable_scip_cuts"] == "true")

    vpc_config["vpolycut_frequency"] = prompt_user(;
        message = "VPolycut Frequency",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        default = "0"
    )
    vpc_config["vpolycut_frequency"] = parse(Int, vpc_config["vpolycut_frequency"])
    vpc_config["vpolycut_priority"] = prompt_user(;
        message = "VPolycut Priority",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        default = "99999"
    )
    vpc_config["vpolycut_priority"] = parse(Int, vpc_config["vpolycut_priority"])
    vpc_config["vpolycut_delay"] = prompt_user(;
        message = "Delay VPolycut (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        default = "false"
    )
    vpc_config["vpolycut_delay"] = (vpc_config["vpolycut_delay"] == "true")

    vpc_config["vpolycut_max_round"] = prompt_user(;
        message = "VPolycut maximum number of cutting plane round to participate in",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        default = "1"
    )
    vpc_config["vpolycut_max_round"] = parse(
        Int, vpc_config["vpolycut_max_round"]
    )

    vpc_config["vpolycut_max_cut_per_round"] = prompt_user(;
        message = "VPolycut Max Cut Per Round",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        default = "150"
    )
    vpc_config["vpolycut_max_cut_per_round"] = parse(
        Int, vpc_config["vpolycut_max_cut_per_round"]
    )

    vpc_config["vpolycut_max_consecutive_fail"] = prompt_user(;
        message = "VPolycut maximum number of fail in cut generation",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        default = "5"
    )
    vpc_config["vpolycut_max_consecutive_fail"] = parse(
        Int, vpc_config["vpolycut_max_consecutive_fail"]
    )

    vpc_config["vpolycut_min_gap_closed_increase"] = prompt_user(;
        message = "VPolycut Minimum Increase Of Disjunctive Gap Closed to be counted as not stagnating",
        validation = (x) ->
            (
                !isnothing(tryparse(Float64, x)) && 0 <= parse(Float64, x) &&
                parse(Float64, x) <= 1
            ),
        error_message = "Not a Float between 0 and 1",
        default = "0.01"
    )
    vpc_config["vpolycut_min_gap_closed_increase"] = parse(
        Float64, vpc_config["vpolycut_min_gap_closed_increase"]
    )

    vpc_config["vpolycut_max_consecutive_stagnation"] = prompt_user(;
        message = "VPolycut maximum number of consecutive cut generations without improvement",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        default = "10"
    )
    vpc_config["vpolycut_max_consecutive_stagnation"] = parse(
        Int, vpc_config["vpolycut_max_consecutive_stagnation"]
    )

    config_file = joinpath(experiment_path, "vpc_config.toml")
    open(config_file, "w") do io
        TOML.print(io, vpc_config)
    end

    # Read the bash template file and write to the experiment folder
    template_file = joinpath(@__DIR__, "job_templates", "vpc_template.sh")
    data = Dict(
        "N" => length(instances),
        "EXPERIMENT_PATH" => experiment_path,
        "CONFIG_PATH" => config_file
    )
    bash_template = read(template_file, String)
    bash_script = Mustache.render(bash_template, data)
    output_file = joinpath(experiment_path, "job_script.sh")
    open(output_file, "w") do f
        write(f, bash_script)
    end
    run(`sbatch $output_file`)
end