using SCIP
using SparseArrays

struct Point <: AbstractVector{SCIP.SCIP_Real}
    coordinates::SparseVector{SCIP.SCIP_Real}
end

function Point(dimension::Int)
    return Point(spzeros(dimension))
end

function Point(coordinates::Vector{SCIP.SCIP_Real})
    point = Point(length(coordinates))
    for i in 1:length(coordinates)
        if !is_zero(coordinates[i])
            point[i] = coordinates[i]
        end
    end
    return point
end

function SparseArrays.findnz(point::Point)
    return findnz(point.coordinates)
end

function Base.getindex(point::Point, i::Int)
    return point.coordinates[i]
end

function Base.setindex!(point::Point, value::SCIP.SCIP_Real, i::Int)
    point.coordinates[i] = value
end

function Base.size(point::Point)
    return size(point.coordinates)
end

function Base.:-(
    point1::Point,
    point2::Point
)::Point
    if length(point1) != length(point2)
        error("Points must have the same dimension")
    end
    temp = Point(point1.coordinates - point2.coordinates)
    # New zero entry may occur so we need to clean the point
    return clean(temp)
end

function clean(point::Point)::Point
    index, _ = findnz(point.coordinates)
    for i in index
        if is_zero(point.coordinates[i])
            point.coordinates[i] = 0
        end
    end
    return point
end

function Base.iterate(point::Point, i::Int = 1)
    if i <= length(point)
        return (point[i], i + 1)
    else
        return nothing
    end
end

function as_dense_vector(point::Point)::Vector{SCIP.SCIP_Real}
    return Vector(point.coordinates)
end