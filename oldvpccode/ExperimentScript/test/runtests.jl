using VPolyCut
using ExperimentScript
using Test
using SCIP
@testset "Basic Test" begin
    cd("..")
    execution_parameters = ExperimentScript.ExecutionParameters(
        "neos5",
        "vpc",
        VPolyCut.IntersectionSeparator,
        false)
    experiment_store = ExperimentScript.setup_environment(execution_parameters)
    ExperimentScript.load_problem_data(experiment_store)
    SCIP.@SCIP_CALL SCIP.SCIPsolve(experiment_store.scip)
    ExperimentScript.record_result(experiment_store)
    @test true
end