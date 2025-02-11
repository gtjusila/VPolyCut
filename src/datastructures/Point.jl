using SCIP

struct Point <: AbstractVector{SCIP.SCIP_Real}
    coordinates::Vector{SCIP.SCIP_Real}
end

function Point(dimension::Int)
    return Point(zeros(dimension))
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
    return Point(point1.coordinates - point2.coordinates)
end

function Base.iterate(point::Point, i::Int = 1)
    if i <= length(point)
        return (point[i], i + 1)
    else
        return nothing
    end
end

function as_vector(point::Point)::Vector{SCIP.SCIP_Real}
    return point.coordinates
end