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
        heuristics = false,
    )
    run_experiment(instance, config)

end

run_test_miplib("test_box")