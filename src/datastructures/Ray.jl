using SCIP
using Lazy

mutable struct Ray <: AbstractVector{SCIP.SCIP_Real}
    coefficients::Vector{SCIP.SCIP_Real}
    generating_variable::Variable
end
function get_coefficients(ray::Ray)
    return ray.coefficients
end

function get_generating_variable(ray::Ray)
    return ray.generating_variable
end

function set_coefficients!(ray::Ray, coefficients::Vector{SCIP.SCIP_Real})
    ray.coefficients = coefficients
end

@forward Ray.coefficients Base.size, Base.getindex, Base.setindex!
