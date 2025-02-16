using SCIP
using Lazy
using SparseArrays

struct Ray <: AbstractVector{SCIP.SCIP_Real}
    coefficients::SparseVector{SCIP.SCIP_Real}
    objective_coefficient::SCIP.SCIP_Real
end

function Ray(dimension::Int, objective_coefficient::SCIP.SCIP_Real)::Ray
    return Ray(spzeros(dimension), objective_coefficient)
end

function SparseArrays.nnz(ray::Ray)
    return nnz(ray.coefficients)
end

function as_dense_vector(ray::Ray)::Vector{SCIP.SCIP_Real}
    return Vector(ray.coefficients)
end

function SparseArrays.findnz(ray::Ray)
    return findnz(ray.coefficients)
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

function Base.size(ray::Ray)
    return size(ray.coefficients)
end