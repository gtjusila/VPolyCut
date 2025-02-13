# Helper structure for managing pointray collection including for removing duplicate rays
using SCIP
@kwdef mutable struct PointRayCollection
    points::Vector{Point} = []
    rays::Vector{Ray} = []
end

get_points(collection::PointRayCollection) = collection.points
get_rays(collection::PointRayCollection) = collection.rays
num_points(collection::PointRayCollection) = length(collection.points)
num_rays(collection::PointRayCollection) = length(collection.rays)

function add_point(
    collection::PointRayCollection,
    new_point::Point
)
    push!(collection.points, new_point)
end

function add_ray(collection::PointRayCollection, new_ray::Ray)
    push!(collection.rays, new_ray)
end

function dimension(collection::PointRayCollection)
    if isempty(collection.points)
        throw("Cannot get dimension of an empty PointRayCollection")
    else
        return length(collection.points[1])
    end
end