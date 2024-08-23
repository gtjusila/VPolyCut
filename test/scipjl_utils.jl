using VPolyhedralCut.SCIPJLUtils
using Test
using JuMP
using SCIP

@testset "Simple Test" begin
    model = setup_scip_safe_jump_model()
    set_mode_cutting_plane_experiment(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    data = Dict()
    data["called"] = 0
    rootlp_callback = initiate_callback(model, data) do data
        data["called"] += 1
    end

    @variable(model, x >= 0, Bin)
    @variable(model, y >= 0, Bin)
    @constraint(model, c1, x + y <= 1.5)
    @objective(model, Min, -x)

    scip = get_scip_data_from_model(model)
    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(scip)
    register_callback(model, SCIP.SCIP_EVENTTYPE_NODESOLVED, rootlp_callback)

    JuMP.optimize!(model)
    @test data["called"] == 1
end