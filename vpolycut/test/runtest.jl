@testitem "Simplex" begin
    import JuMP
    include("utilities.jl")
    model = setup_jump_model()
    scip = get_scipdata_from_model(model)
    include_separator(scip, VPolyCut.IntersectionSeparator; scipd=scip)
    set_scip_parameters_easy(model)

    JuMP.@variable(model, x >= 0, Int)
    JuMP.@variable(model, y >= 0, Int)
    JuMP.@constraint(model, x + y <= 2.5)
    JuMP.@objective(model, Max, 2 * x + y)

    JuMP.optimize!(model)
end