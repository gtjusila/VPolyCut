using SCIP
using VPolyhedralCut
using VPolyhedralCut.SCIPJLUtils
using JuMP
using Test

# Test Branch and Bound
# https://www.ie.bilkent.edu.tr/~mustafap/courses/bb.pdf
@testset "Branch and bound test 1" begin
    model = setup_scip_safe_jump_model()
    set_everything_off(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    @variable(model, 5 >= x >= 0, Int)
    @variable(model, y >= 0, Int)
    @constraint(model, c1, -10 * x + 20 * y <= 22)
    @constraint(model, c2, 5 * x + 10 * y <= 49)
    @objective(model, Max, -x + 4 * y)

    data = Dict()
    data["scip"] = get_scip_data_from_model(model)
    data["executed"] = Ref{Bool}(false)

    firstlp = initiate_callback(model, data) do data
        scip = data["scip"]
        branchandbound = VPolyhedralCut.BranchAndBound(scip)

        # ILP is feasible
        @test VPolyhedralCut.execute_branchandbound(branchandbound) == true
        x = VPolyhedralCut.get_best_solution(branchandbound)
        # Check Solution
        @test x[1] == 2 && x[2] == 2

        # Obj Sense is Max so reverse
        obj = VPolyhedralCut.get_primal_bound(branchandbound)
        @test obj == -6

        SCIP.SCIPinterruptSolve(scip)
        data["executed"][] = true
    end

    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(data["scip"])
    register_callback(model, SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, firstlp)

    optimize!(model)
    @test data["executed"][]
end

# Test Branch and Bound
# https://www.ie.bilkent.edu.tr/~mustafap/courses/bb.pdf
@testset "Branch and bound test 2" begin
    model = setup_scip_safe_jump_model()
    set_everything_off(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    @variable(model, x[1:4], Bin)
    @constraint(model, c1, 6 * x[1] + 3 * x[2] + 5 * x[3] + 2 * x[4] <= 10)
    @constraint(model, c2, x[3] + x[4] <= 1)
    @constraint(model, c3, -x[1] + x[3] <= 0)
    @constraint(model, c4, -x[2] + x[4] <= 0)
    @objective(model, Max, 9 * x[1] + 5 * x[2] + 6 * x[3] + 4 * x[4])

    data = Dict()
    data["scip"] = get_scip_data_from_model(model)
    data["executed"] = Ref{Bool}(false)
    firstlp = initiate_callback(model, data) do data
        scip = data["scip"]
        branchandbound = VPolyhedralCut.BranchAndBound(scip)

        # ILP is feasible
        @test VPolyhedralCut.execute_branchandbound(branchandbound) == true
        x = VPolyhedralCut.get_best_solution(branchandbound)
        # Check Solution
        @test x[1] == 1 && x[2] == 1 && x[3] == 0 && x[4] == 0

        # Obj Sense is Max so reverse
        obj = VPolyhedralCut.get_primal_bound(branchandbound)
        @test obj == -14

        SCIP.SCIPinterruptSolve(scip)
        data["executed"][] = true
    end

    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(data["scip"])
    register_callback(model, SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, firstlp)

    optimize!(model)
    @test data["executed"][]
end

# http://web.tecnico.ulisboa.pt/mcasquilho/compute/_linpro/TaylorB_module_c.pdf
@testset "Branch and bound test 3" begin
    model = setup_scip_safe_jump_model()
    set_everything_off(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    @variable(model, x >= 0, Int)
    @variable(model, y >= 0, Int)
    @constraint(model, c1, 2 * x + y <= 10)
    @constraint(model, c2, 3 * x + 6 * y <= 40)
    @objective(model, Max, 2 * x + 3 * y)

    data = Dict()
    data["scip"] = get_scip_data_from_model(model)
    data["executed"] = Ref{Bool}(false)

    firstlp = initiate_callback(model, data) do data
        scip = data["scip"]
        branchandbound = VPolyhedralCut.BranchAndBound(scip)

        # ILP is feasible
        @test VPolyhedralCut.execute_branchandbound(branchandbound) == true

        x = VPolyhedralCut.get_best_solution(branchandbound)
        # Check Solution
        @test x[1] == 1 && x[2] == 6

        # Obj Sense is Max so reverse
        obj = VPolyhedralCut.get_primal_bound(branchandbound)
        @test obj == -20

        SCIP.SCIPinterruptSolve(scip)
        data["executed"][] = true
    end

    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(data["scip"])
    register_callback(model, SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, firstlp)

    optimize!(model)
    @test data["executed"][]
end

# https://ocw.mit.edu/courses/15-053-optimization-methods-in-management-science-spring-2013/1f26dafa0483e3283ff5a3774f7fc194_MIT15_053S13_lec12.pdf
@testset "Branch and bound test 4" begin
    model = setup_scip_safe_jump_model()
    set_everything_off(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    @variable(model, x[1:4] >= 0, Bin)
    @constraint(model, c1, 8 * x[1] + x[2] + 5 * x[3] + 4 * x[4] <= 9)
    @objective(model, Max, 12 * x[1] + x[2] + 10 * x[3] + 2 * x[4])

    data = Dict()
    data["scip"] = get_scip_data_from_model(model)
    data["executed"] = Ref{Bool}(false)

    firstlp = initiate_callback(model, data) do data
        scip = data["scip"]
        branchandbound = VPolyhedralCut.BranchAndBound(scip)

        # Check that is feasible
        @test VPolyhedralCut.execute_branchandbound(branchandbound) == true

        x = VPolyhedralCut.get_best_solution(branchandbound)
        # Check Solution
        @test x[1] == 1 && x[2] == 1 && x[3] == 0 && x[4] == 0

        # Obj Sense is Max so reverse
        obj = VPolyhedralCut.get_primal_bound(branchandbound)
        @test obj == -13

        SCIP.SCIPinterruptSolve(scip)
        data["executed"][] = true
    end

    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(data["scip"])
    register_callback(model, SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, firstlp)

    optimize!(model)
    @test data["executed"][]
end