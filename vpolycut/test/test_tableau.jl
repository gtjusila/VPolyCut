@testitem "Tableau_Simplex" begin
    import JuMP
    import SCIP
    include("utilities.jl")
    include("test_tableau_helper.jl")

    model = setup_jump_model()
    scip = get_scipdata_from_model(model)
    include_separator(scip, LPTableau; scipd=scip)
    set_scip_parameters_easy(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    JuMP.@variable(model, x >= 0, Int)
    JuMP.@variable(model, y >= 0, Int)
    JuMP.@constraint(model, x + y <= 2.5)
    JuMP.@objective(model, Max, 2 * x + y)

    JuMP.optimize!(model)
end

@testitem "Tableau_TestCase_1" begin
    using JuMP
    import SCIP
    include("utilities.jl")
    include("test_tableau_helper.jl")

    model = setup_jump_model()
    scip = get_scipdata_from_model(model)
    include_separator(scip, LPTableau; scipd=scip)
    set_scip_parameters_easy(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    JuMP.@variable(model, x1 >= 0, Int)
    JuMP.@variable(model, x2 >= 0, Int)
    JuMP.@variable(model, x3 >= 0, Int)
    JuMP.@variable(model, x4 >= 0, Int)

    # Define the constraints
    @constraint(model, 2x1 + x2 + x3 + x4 <= 7)
    @constraint(model, x1 + 3x2 + 2x3 + 3x4 <= 12)
    @constraint(model, 2x1 + 2x2 + 3x3 + x4 <= 10)


    JuMP.@objective(model, Max, 2 * x1 + 3 * x2 + 4 * x3 + x4)

    print(model)
    JuMP.optimize!(model)
end

@testitem "Tableau_Conforti" begin
    using JuMP
    import SCIP
    include("utilities.jl")
    include("test_tableau_helper.jl")

    model = setup_jump_model()
    scip = get_scipdata_from_model(model)
    include_separator(scip, LPTableau; scipd=scip)
    set_scip_parameters_easy(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    JuMP.@variable(model, x1 >= 0, Int)
    JuMP.@variable(model, x2 >= 0, Int)
    JuMP.@variable(model, x3 >= 0, Int)

    # Define the objective function
    @objective(model, Min, -0.5x2 + -x3)

    # Define the constraints
    @constraint(model, c1, x1 + x2 + x3 <= 2)
    @constraint(model, c2, x1 - 0.5 * x3 >= 0)
    @constraint(model, c3, x2 - 0.5 * x3 >= 0)
    @constraint(model, c4, x1 + 0.5 * x3 <= 1)
    @constraint(model, c5, -x1 + x2 + x3 <= 1)

    print(model)
    JuMP.optimize!(model)
end