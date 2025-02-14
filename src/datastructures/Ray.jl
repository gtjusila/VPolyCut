using SCIP
using Lazy
using SparseArrays

struct Ray <: AbstractVector{SCIP.SCIP_Real}
    coefficients::SparseVector{SCIP.SCIP_Real}
    generating_variable::Union{Variable,Nothing}
end

function Ray(dimension::Int, generating_variable::Union{Variable,Nothing})
    return Ray(spzeros(dimension), generating_variable)
end

function Ray(dimension::Int)
    return Ray(dimension, nothing)
end

function Ray(coefficients::Vector{SCIP.SCIP_Real}, generating_variable::Variable)
    ray = Ray(length(coefficients), generating_variable)
    for i in 1:length(coefficients)
        if !is_zero(coefficients[i])
            ray[i] = coefficients[i]
        end
    end
    return ray
end

function SparseArrays.nnz(ray::Ray)
    return nnz(ray.coefficients)
end

function Ray(coefficients::Vector{SCIP.SCIP_Real})
    return Ray(coefficients, Variable())
end

function as_dense_vector(ray::Ray)
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