using DataFrames
using CSV
using Mustache
include("utils.jl")

println("Setup experiment")

instance_list = prompt_user(;
    message="Instance List",
    validation=(x) -> isfile(abspath(x)),
    error_message="Invalid Path.",
    default="experiment_data/instances/instances_list.txt"
)

instance_dir = prompt_user(;
    message="Instance List",
    validation=(x) -> isdir(abspath(x)),
    error_message="Invalid Path.",
    default="experiment_data/instances"
)

mode = prompt_user(;
    message="Mode",
    validation=(x) -> x in ["vpc", "gomory"],
    error_message="Invalid Mode.",
    default="vpc"
)

# If experiment_runs is not there create it first.
runs_path = abspath(joinpath(dirname(@__FILE__), "..", "experiment_runs"))
if !isdir(runs_path)
    @info "Creating experiment_runs folder"
    mkdir(runs_path)
end

experiment_path = prompt_user(;
    message="Experiment Name:",
    validation=(x) ->
        ((x != "") && !isdir(joinpath(runs_path, x))),
    error_message="Invalid Mode.",
    default=""
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
        id=1:length(instances),
        instances_path=instances_path,
        output_path=output_path
    )
    run_settings_file = joinpath(experiment_path, "experiment_list.tsv")
    CSV.write(run_settings_file, config_dataframe; delim='\t')

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
        id=1:length(instances),
        instances_path=instances_path,
        output_path=output_path
    )
    run_settings_file = joinpath(experiment_path, "experiment_list.tsv")
    CSV.write(run_settings_file, config_dataframe; delim='\t')

    for path in instances_path
        if (!isfile(path))
            @warn "Possibly missing instance $(path)"
        end
    end

    # Read the bash template file and write to the experiment folder
    template_file = joinpath(@__DIR__, "job_templates", "vpc_template.sh")
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
end