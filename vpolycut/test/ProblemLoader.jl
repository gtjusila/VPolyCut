using JuMP
abstract type ProblemLoader end

function load_problem(loader::T, model) where {T<:ProblemLoader}
    error("load_problem not implemented for type ", T)
end

struct ProblemConforti <: ProblemLoader end

function load_problem(loader::ProblemConforti, model)
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
end