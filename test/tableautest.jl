using SCIP
using VPolyhedralCut
using VPolyhedralCut.SCIPJLUtils
using JuMP

# Test the Retrieval of various tableau information
@testset "Test Tableau Info" begin
    model = setup_scip_safe_jump_model()
    set_mode_cutting_plane_experiment(model)
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
            tableau = VPolyhedralCut.construct_tableau(scip)
            col = VPolyhedralCut.get_var_from_column(tableau, 1)
            @test VPolyhedralCut.get_lb(col) == 0.0
            @test VPolyhedralCut.get_ub(col) == SCIP.SCIPinfinity(scip)
            @test VPolyhedralCut.get_basis_status(col) == SCIP.SCIP_BASESTAT_BASIC
            @test VPolyhedralCut.get_noriginalcols(tableau) == 2
            @test VPolyhedralCut.get_noriginalrows(tableau) == 1
            @test VPolyhedralCut.get_nvars(tableau) == 3
            @test VPolyhedralCut.get_nbasis(tableau) == 1
            @test VPolyhedralCut.get_sol(col) == 1.5
            row = VPolyhedralCut.get_var_from_column(tableau, 3)
            @test VPolyhedralCut.get_lb(row) == -1.5
            @test VPolyhedralCut.get_basis_status(row) == SCIP.SCIP_BASESTAT_LOWER
            @test VPolyhedralCut.get_ub(row) == SCIP.SCIPinfinity(scip)
            @test VPolyhedralCut.get_sol(row) == -1.5
        end
    end

    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(scipd)
    register_callback(model, SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, firstlp_callback)

    JuMP.optimize!(model)
end
