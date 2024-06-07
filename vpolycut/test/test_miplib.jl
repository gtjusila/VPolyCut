using Test
import Polyhedra
import SCIP
import MathOptInterface as MOI

include("./helper_experiment.jl")
include("./helper_intersection.jl")

function run_test_miplib(instance::String) 

    config = ExperimentConfiguration(
        debug = true,
        separator = false,
        heuristics = true,
        vpolycut = false,
        gomory =  true,
        vpolycut_limit = 20,
        node_limit=1 
    )
    run_experiment(instance, config)

end

run_test_miplib("gen-ip054")