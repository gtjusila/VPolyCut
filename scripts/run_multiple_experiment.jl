using ArgParse
using TOML
using Printf
using Dates
using DataFrames
using CSV
using Mustache

# Get Global Config
global_config = open(joinpath(@__DIR__, "global_config.toml")) do file
    return TOML.parse(file)
end
instances_base_path = global_config["instances_base_path"]

# Setup and get CLI arguments
settings = ArgParseSettings()
@add_arg_table settings begin
    "--config_path", "-c"
    help = "path to TOML configuration file for the experiment"
    required = true
end
cli_args = parse_args(settings)
config_path = abspath(cli_args["config_path"])

# Parse the configuration file
experiment_config = open(config_path) do file
    return TOML.parse(file)
end
if isnothing(experiment_config)
    error("Could not parse the configuration file")
end

# Prepare Experiment Runs Folder
runs_path = joinpath(pwd(), "experiment_runs")
if !isdir(runs_path)
    mkdir(runs_path)
end

# Prepare the run folder for this experiment
now = Dates.now()
timestamp = Dates.format(now, "yyyymmddTHHMMSS")
folder_name = string(timestamp)
experiment_path = joinpath(runs_path, timestamp)
@info "Experiment folder will be done at $experiment_path"
mkdir(experiment_path)

# Prepare the instances  
modes = experiment_config["modes"]
instances = experiment_config["instances"]

# Generate all combinations of modes and instances
mode_instance_combinations = vec(collect(Iterators.product(modes, instances)))
label = label = [item[2] * "_" * item[1] for item in mode_instance_combinations]
mode = [item[1] for item in mode_instance_combinations]
instance_path = [
    instances_base_path * item[2] * ".mps" for item in mode_instance_combinations
]
output_path = [experiment_path * "/" * item for item in label]
for path in output_path
    if !isdir(path)
        mkdir(path)
    end
end

# Create a DataFrame with ID, mode, and instance
config_dataframe = DataFrame(;
    id=1:length(mode_instance_combinations),
    label=label,
    mode=mode,
    instance_path=instance_path,
    output_path=output_path
)

# Write the DataFrame to a tab-delimited text file
run_settings_file = joinpath(experiment_path, "experiment_list.tsv")
CSV.write(run_settings_file, config_dataframe; delim='\t')

# Read the bash template file and write to the experiment folder
template_file = joinpath(@__DIR__, "job_template.sh")
data = Dict(
    "N" => length(mode_instance_combinations),
    "JULIA_DEPOT_PATH" => global_config["julia_depot_path"],
    "PATH_TO_SCRIPT" => global_config["path_to_script"]
)
bash_template = read(template_file, String)
bash_script = Mustache.render(bash_template, data)
output_file = joinpath(experiment_path, "job_script.sh")
open(output_file, "w") do f
    write(f, bash_script)
end

runs_path = joinpath(pwd(), "experiment_runs")
if !isdir(runs_path)
    mkdir(runs_path)
end
