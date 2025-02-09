using SCIP
using Lazy

mutable struct Ray <: AbstractVector{SCIP.SCIP_Real}
    coefficients::Vector{SCIP.SCIP_Real}
    generating_variable::Variable
end
function Ray()
    return Ray([], Variable())
end
function Ray(coefficients::Vector{SCIP.SCIP_Real})
    return Ray(coefficients, Variable())
end
function get_coefficients(ray::Ray)
    return ray.coefficients
end

function Base.getindex(ray::Ray, i::Int)
    return ray.coefficients[i]
end

function get_generating_variable(ray::Ray)
    return ray.generating_variable
end

function Base.setindex!(ray::Ray, value::SCIP.SCIP_Real, i::Int)
    ray.coefficients[i] = value
end

function set_coefficients!(ray::Ray, coefficients::Vector{SCIP.SCIP_Real})
    ray.coefficients = coefficients
end

@forward Ray.coefficients Base.size, Base.getindex, Base.setindex!
