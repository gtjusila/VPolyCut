using Test
import Polyhedra
import SCIP
import MathOptInterface as MOI

include("./experiment_runexperiment.jl")
include("./helper_intersection.jl")

"""
Setup experiment folder. The format is xxxxx_instancename_date where xxxxx is a random
alphanumeric code
"""
function setup_experiment_folder(instance::String)
    if !isdir("experiments")
        mkdir("experiments")
    end
    code = lowercase(randstring(5))
    date = Dates.format(Dates.today(), "yyyymmdd")
    path = joinpath(pwd(), "experiments", "$(code)_$(instance)_$(date)")
    mkdir(path)
    string = "Created Experiment Folder $(path)"
    seperator = join(["=" for i=1:length(string)+10])
    println(seperator)
    println(string)
    println(seperator)
    return path
end

function run_test_miplib(instance::String) 
    path = setup_experiment_folder(instance)
    config = ExperimentConfiguration(
        path = path,
        type = "vpoly", 
        debug = true,
        separator = false,
        heuristics = false,
        vpolycut = false,
        gomory =  true,
        vpolycut_limit = 33,
        node_limit=1,
        max_separounds=1
    )
    run_experiment(instance, config)
end

run_test_miplib("pg")