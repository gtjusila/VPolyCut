include("utilities.jl")
include("utilities/tableau_test_utilities.jl")
include("ProblemLoader.jl")
using Revise
using JuMP
using VPolyCut
using Test

@testset "Tableau Test" begin
    model = setup_jump_model()
    scip = get_scipdata_from_model(model)
    function testcase(scip)
        tableau = VPolyCut.construct_tableau(scip)
        col = VPolyCut.get_var_from_column(tableau, 1)
        @test VPolyCut.get_lb(col) == 0.0
        @test VPolyCut.get_ub(col) == SCIP.SCIPinfinity(scip)
        @test VPolyCut.get_noriginalcols(tableau) == 3
        @test VPolyCut.get_noriginalrows(tableau) == 5
        @test VPolyCut.get_nvars(tableau) == 8
        @test VPolyCut.get_nbasis(tableau) == 5
        row = VPolyCut.get_var_from_column(tableau, 4)
        @test VPolyCut.get_lb(row) == -SCIP.SCIPinfinity(scip)
        @test VPolyCut.get_ub(row) == 2
    end
    sepa = include_separator(scip, TableauTest; scipd=scip, callback=testcase)
    set_scip_parameters_easy(model)
    problem = ProblemConforti()
    load_problem(problem, model)
    JuMP.optimize!(model)
end