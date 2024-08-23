@testitem "Simplex" begin
    import JuMP
    include("../utilities.jl")
    model = setup_jump_model()
    scip = get_scipdata_from_model(model)
    sepa = include_separator(scip, VPolyCut.IntersectionSeparator; scipd=scip)
    set_scip_parameters_easy(model)

    JuMP.@variable(model, x >= 0, Int)
    JuMP.@variable(model, y >= 0, Int)
    JuMP.@constraint(model, x + y <= 2.5)
    JuMP.@objective(model, Max, 2 * x + y)

    JuMP.optimize!(model)
end

@testitem "Trapezoid" default_imports = false begin
    using Revise
    using JuMP
    using VPolyCut
    include("../utilities.jl")
    model = setup_jump_model()
    scip = get_scipdata_from_model(model)
    sepa = include_separator(scip, VPolyCut.IntersectionSeparator; scipd=scip)
    set_scip_parameters_easy(model)

    JuMP.@variable(model, x1 >= 0, Int)
    JuMP.@variable(model, x2 >= 0, Int)

    # Define the objective function
    @objective(model, Max, x1 + x2)

    # Define the constraints
    @constraint(model, c1, x1 + 0.00000001 * x2 <= 2.25)
    @constraint(model, c2, x1 - 10 * x2 >= -25)

    JuMP.optimize!(model)
end

@testitem "Conforti" default_imports = false begin
    using Revise
    using JuMP
    using VPolyCut
    include("../utilities.jl")
    model = setup_jump_model()
    scip = get_scipdata_from_model(model)
    sepa = include_separator(scip, VPolyCut.IntersectionSeparator; scipd=scip)
    set_scip_parameters_easy(model)

    JuMP.@variable(model, x1 >= 0, Int)
    JuMP.@variable(model, x2 >= 0, Int)
    JuMP.@variable(model, x3 >= 0, Int)

    # Define the objective function
    @objective(model, Max, 0.5x2 + x3)

    # Define the constraints
    @constraint(model, c1, x1 + x2 + x3 <= 2)
    @constraint(model, c2, x1 - 0.5 * x3 >= 0)
    @constraint(model, c3, x2 - 0.5 * x3 >= 0)
    @constraint(model, c4, x1 + 0.5 * x3 <= 1)
    @constraint(model, c5, -x1 + x2 + x3 <= 1)

    JuMP.optimize!(model)
end