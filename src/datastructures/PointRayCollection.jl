# Helper structure for managing pointray collection including for removing duplicate rays
using SCIP

"""
    CornerPoint

A CornerPoint represents a vertex of the LP. It have a point which holds the actual coordinates 
and objective value associated to it 
"""
@kwdef struct CornerPoint
    point::Point
    objective_value::SCIP.SCIP_Real
    orig_objective_value::SCIP.SCIP_Real
end

function get_point(corner::CornerPoint)::Point
    return corner.point
end

function get_objective_value(corner::CornerPoint)::SCIP.SCIP_Real
    return corner.objective_value
end

function get_orig_objective_value(corner::CornerPoint)::SCIP.SCIP_Real
    return corner.orig_objective_value
end

@forward CornerPoint.point Base.size, Base.getindex, Base.setindex!

@kwdef mutable struct PointRayCollection
    points::Vector{CornerPoint} = []
    rays::Vector{Ray} = []
end

get_points(collection::PointRayCollection) = collection.points
get_rays(collection::PointRayCollection) = collection.rays
num_points(collection::PointRayCollection) = length(collection.points)
num_rays(collection::PointRayCollection) = length(collection.rays)

function add_point(
    collection::PointRayCollection,
    new_point::Point,
    objective_value::SCIP.SCIP_Real,
    orig_objective_value::SCIP.SCIP_Real
)
    push!(collection.points, CornerPoint(new_point, objective_value, orig_objective_value))
end

function add_ray(collection::PointRayCollection, new_ray::Ray)
    push!(collection.rays, new_ray)
end

"""
A helper function to remove duplicate ray
"""
function remove_duplicate_rays(collection::PointRayCollection)
    rays = get_rays(collection)
    unique_indices = get_unique_row_indices(
        map(ray -> get_coefficients(ray), collection.rays),
        collection.scip
    )
    collection.rays = rays[unique_indices]
end