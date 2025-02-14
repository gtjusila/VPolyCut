using SCIP
using SparseArrays

struct Ray <: AbstractVector{SCIP.SCIP_Real}
    coefficients::SparseVector{SCIP.SCIP_Real}
    objective_coefficient::SCIP.SCIP_Real
end

function Ray(dimension::Int, objective_coefficient::SCIP.SCIP_Real)
    return Ray(spzeros(dimension), objective_coefficient)
end

function get_objective(ray::Ray)
    return ray.objective_coefficient
end

function as_dense_vector(ray::Ray)
    return Vector(ray.coefficients)
end

function SparseArrays.nnz(ray::Ray)
    return nnz(ray.coefficients)
end

function SparseArrays.findnz(ray::Ray)
    return findnz(ray.coefficients)
end

function Base.getindex(ray::Ray, i::Int)
    return ray.coefficients[i]
end

function Base.setindex!(ray::Ray, value::SCIP.SCIP_Real, i::Int)
    ray.coefficients[i] = value
end

function Base.size(ray::Ray)
    return size(ray.coefficients)
end