using SCIP
using SparseArrays

struct Point <: AbstractVector{SCIP.SCIP_Real}
    coordinates::SparseVector{SCIP.SCIP_Real}
    orig_objective_value::SCIP.SCIP_Real
end

### Constructors ###
function Point(dimension::Int, orig_objective_value::SCIP.SCIP_Real)
    return Point(spzeros(dimension), orig_objective_value)
end

function Point(coordinates::Vector{SCIP.SCIP_Real})
    return Point(coordinates, 0.0)
end

### Methods ###

function clean(point::Point)::Point
    index, _ = findnz(point.coordinates)
    for i in index
        if is_zero(point.coordinates[i])
            point.coordinates[i] = 0
        end
    end
    return point
end

function as_dense_vector(point::Point)::Vector{SCIP.SCIP_Real}
    return Vector(point.coordinates)
end

function get_objective_value(point::Point)::SCIP.SCIP_Real
    return point.orig_objective_value
end

### Temporary Solutions ###
function set_objective_value!(point::Point, value::SCIP.SCIP_Real)::Point
    return Point(point.coordinates, value)
end
function substract(point1::Point, point2)::Point
    temp = Point(
        point1.coordinates - point2.coordinates,
        point1.orig_objective_value
    )
    clean(temp)
    return temp
end

### Interface to SparseArrays ###

function SparseArrays.findnz(point::Point)
    return findnz(point.coordinates)
end

function SparseArrays.nnz(point::Point)
    return nnz(point.coordinates)
end

### Interface to Base.AbstractVector ###

function Base.getindex(point::Point, i::Int)
    return point.coordinates[i]
end

function Base.setindex!(point::Point, value::SCIP.SCIP_Real, i::Int)
    point.coordinates[i] = value
end

function Base.size(point::Point)
    return size(point.coordinates)
end

function Base.iterate(point::Point, i::Int = 1)
    if i <= length(point)
        return (point[i], i + 1)
    else
        return nothing
    end
end
