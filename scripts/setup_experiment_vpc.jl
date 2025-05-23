using DataFrames
using CSV
using Mustache
using TOML

include("utils.jl")
println("----------------------")
println("Setup VPC experiment")
println("----------------------")
println("")
println("------------------")
println("Instance Setting")
println("------------------")
instance_list = prompt_user(;
    message = "Instance List",
    validation = (x) -> isfile(abspath(x)),
    error_message = "Invalid Path.",
    default = "experiment_data/miplibbench/instances_list.txt"
)

instance_dir = prompt_user(;
    message = "Instance Dir",
    validation = (x) -> isdir(abspath(x)),
    error_message = "Invalid Path.",
    default = "experiment_data/miplibbench"
)

solution_dir = prompt_user(;
    message = "Solution Directory",
    validation = (x) -> isdir(abspath(x)),
    error_message = "Invalid Path.",
    default = "experiment_data/miplibbenchsolutions"
)

use_solution = prompt_user(;
    message = "Use solution",
    validation = (x) -> x == "true" || x == "false",
    error_message = "Neither true nor false.",
    parse = (x) -> x == "true",
    default = "true"
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
    parse = (x) -> joinpath(runs_path, x),
    default = ""
)
mkdir(experiment_path)

### Done with basic parameter, start reading input ###
instances = readlines(instance_list)

### Creating Config TSV ###
instances_path = [
    abspath(joinpath(instance_dir, instance * ".mps")) for instance in instances
]
output_path = [abspath(joinpath(experiment_path, instance)) for instance in instances]
solution_path = map(instances) do instance
    if isfile(joinpath(solution_dir, instance * ".sol")) && use_solution
        return abspath(joinpath(solution_dir, instance * ".sol"))
    else
        return ""
    end
end
config_dataframe = DataFrame(;
    id = 1:length(instances),
    instances_path = instances_path,
    output_path = output_path,
    solution_path = solution_path
)
run_settings_file = joinpath(experiment_path, "experiment_list.tsv")
CSV.write(run_settings_file, config_dataframe; delim = '\t')

# Make Path for all instances
for path in output_path
    mkdir(path)
end

# Check if instance are all there
for path in instances_path
    if (!isfile(path))
        @warn "Possibly missing instance $(path)"
    end
end

use_existing_config = prompt_user(;
    message = "Use Existing Config File?",
    validation = (x) -> x == "true" || x == "false",
    error_message = "Neither true nor false.",
    parse = (x) -> x == "true",
    default = "true"
)
if !use_existing_config
    # Create config file
    vpc_config = Dict()
    println()
    println("------------------")
    println("SCIP Settings")
    println("------------------")
    vpc_config["scip_node_limit"] = prompt_user(;
        message = "Node Limit",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "-1"
    )
    vpc_config["scip_max_root_cutting_plane_rounds"] = prompt_user(;
        message = "Max Cutting Plane Rounds on the Root Node",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "-1"
    )
    vpc_config["scip_time_limit"] = prompt_user(;
        message = "Time Limit",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = string(1209600)
    )

    vpc_config["scip_disable_scip_cuts"] = prompt_user(;
        message = "Disable SCIP Cuts (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "true"
    )

    vpc_config["scip_enable_cut_selection"] = prompt_user(;
        message = "SCIP enable cut selection (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "true"
    )

    vpc_config["scip_enable_strong_branching_lookahead"] = prompt_user(;
        message = "SCIP enable strong branching lookahead (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "true"
    )

    vpc_config["scip_enable_heuristic"] = prompt_user(;
        message = "SCIP enable heuristic (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "true"
    )

    vpc_config["scip_enable_conflict_analysis"] = prompt_user(;
        message = "SCIP enable conflict analysis (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "true"
    )

    vpc_config["scip_enable_root_node_propagation"] = prompt_user(;
        message = "SCIP enable root node propagation (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "true"
    )

    vpc_config["scip_allow_restart"] = prompt_user(;
        message = "SCIP allow restart (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "true"
    )

    println()
    println("------------------")
    println("VPC Settings")
    println("------------------")
    vpc_config["vpc_frequency"] = prompt_user(;
        message = "VPolycut Frequency",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "0"
    )

    vpc_config["vpc_priority"] = prompt_user(;
        message = "VPolycut Priority",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "99999"
    )

    vpc_config["vpc_delayed"] = prompt_user(;
        message = "Delay VPolycut (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "false"
    )

    vpc_config["vpc_max_participating_round"] = prompt_user(;
        message = "Maximum number of cutting plane round to participate in",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "1"
    )

    vpc_config["vpc_max_cut_per_round"] = prompt_user(;
        message = "Max Cut Per Round",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "150")

    vpc_config["vpc_min_restart"] = prompt_user(;
        message = "Minimum number of restart before VPC participate",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "1"
    )

    vpc_config["vpc_n_leaves"] = prompt_user(;
        message = "Number of Leaves in Branch and Bound Tree",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "64"
    )

    vpc_config["vpc_prlp_solve_method"] = prompt_user(;
        message = "PRLP Solve Method (1: PRIMAL SIMPLEX, 2:DUAL SIMPLEX, 3: BARRIER)",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "1"
    )

    vpc_config["vpc_prlp_allow_warm_start"] = prompt_user(;
        message = "Allow PRLP warm start of basis (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "true"
    )

    vpc_config["vpc_prlp_apply_beta_scaling"] = prompt_user(;
        message = "Use beta scaling in PRLP (true/false)",
        validation = (x) -> x == "true" || x == "false",
        error_message = "Neither true nor false.",
        parse = (x) -> x == "true",
        default = "false"
    )

    vpc_config["vpc_prlp_max_consecutive_fail"] = prompt_user(;
        message = "Maximum number of consecutive fail in cut generation",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "5"
    )

    vpc_config["vpc_prlp_min_gap_closed_increase"] = prompt_user(;
        message = "Minimum Increase Of Disjunctive Gap Closed to be counted as not stagnating",
        validation = (x) ->
            (
                !isnothing(tryparse(Float64, x)) && 0 <= parse(Float64, x) &&
                parse(Float64, x) <= 1
            ),
        error_message = "Not a Float between 0 and 1",
        parse = (x) -> parse(Float64, x),
        default = "0.01"
    )

    vpc_config["vpc_prlp_max_consecutive_stagnation"] = prompt_user(;
        message = "Maximum number of consecutive cut generations without improvement",
        validation = (x) -> !isnothing(tryparse(Int, x)),
        error_message = "Not an integer.",
        parse = (x) -> parse(Int, x),
        default = "10"
    )

    config_file = joinpath(experiment_path, "vpc_config.toml")
    open(config_file, "w") do io
        TOML.print(io, vpc_config)
    end
else
    config_file_path = prompt_user(;
        message = "Config File Path",
        validation = (x) -> isfile(abspath(x)),
        error_message = "Invalid Path.",
        default = "config_files/default_config.toml"
    )
    config_file = joinpath(experiment_path, "vpc_config.toml")
    cp(config_file_path, config_file; force = true)
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

# Run Job
#run(`sbatch $output_file`)
