using DataFrames
using CSV
using Mustache

println("Setup experiment")

instance_list = ""
while true
    print("Instance List Path: ")
    global instance_list = abspath(readline())
    if (isfile(instance_list))
        break
    else
        println("Invalid Path.")
    end
end

instance_dir = ""
while true
    print("Instance Directory: ")
    global instance_dir = abspath(readline())
    if (isdir(instance_dir))
        break
    else
        println("Invalid Path.")
    end
end

mode = ""
while true
    print("Mode: ")
    global mode = readline()
    if mode == "gomory" || mode == "vpc"
        break
    else
        println("Invalid Mode.")
    end
end

runs_path = abspath(joinpath(dirname(@__FILE__), "..", "experiment_runs"))
if !isdir(runs_path)
    @info "Creating experiment_runs folder"
    mkdir(runs_path)
end

experiment_path = ""
while true
    print("Experiment Name: ")
    experiment_name = readline()
    global experiment_path = abspath(joinpath(runs_path, experiment_name))
    if !isdir(experiment_path)
        mkdir(experiment_path)
        break
    end
end

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
    template_file = joinpath(@__DIR__, "vpc_template.sh")
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