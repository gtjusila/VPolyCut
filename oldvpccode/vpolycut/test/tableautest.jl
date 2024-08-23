using SCIP
using VPolyCut
using JuMP
using SCIPExperimentUtils

# Test the Retrieval of various tableau information
@testset "Test Tableau Info" begin
    model = setup_scip_safe_jump_model()
    SCIPExperimentUtils.set_mode_cutting_plane_experiment(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    @variable(model, x >= 0, Int)
    @variable(model, y >= 0, Int)
    @constraint(model, c1, x + y <= 1.5)
    @objective(model, Min, -x)

    data = Dict()
    scipd = get_scip_data_from_model(model)
    data["scip"] = scipd

    firstlp_callback = initiate_callback(model, data) do data
        scip = data["scip"]
        if SCIP.SCIPgetDepth(scip) == 0
            tableau = VPolyCut.construct_tableau(scip)
            col = VPolyCut.get_var_from_column(tableau, 1)
            @test VPolyCut.get_lb(col) == 0.0
            @test VPolyCut.get_ub(col) == SCIP.SCIPinfinity(scip)
            @test VPolyCut.get_basis_status(col) == SCIP.SCIP_BASESTAT_BASIC
            @test VPolyCut.get_noriginalcols(tableau) == 2
            @test VPolyCut.get_noriginalrows(tableau) == 1
            @test VPolyCut.get_nvars(tableau) == 3
            @test VPolyCut.get_nbasis(tableau) == 1
            @test VPolyCut.get_sol(col) == 1.5
            row = VPolyCut.get_var_from_column(tableau, 3)
            @test VPolyCut.get_lb(row) == -1.5
            @test VPolyCut.get_basis_status(row) == SCIP.SCIP_BASESTAT_LOWER
            @test VPolyCut.get_ub(row) == SCIP.SCIPinfinity(scip)
            @test VPolyCut.get_sol(row) == -1.5
        end
    end

    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(scipd)
    register_callback(model, SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, firstlp_callback)

    JuMP.optimize!(model)
end
