using SCIP

"""
    ObjectiveFunction

A datastructure representing an objective direction
"""
struct ObjectiveFunction
    """
    The objective direction
    """
    direction::Vector{SCIP.SCIP_Real}
    """
    Time limit to solve the PRLP
    """
    time_limit::Float64
    """
    label
    """
    label::String
    """
    Optimal Basis can be saved and reused
    """
    savable::Bool
end

# Constructor methods for different objective functions

function feasibility_objective(prlp::PRLP)::ObjectiveFunction
    return ObjectiveFunction(
        zeros(SCIP.SCIP_Real, prlp.dimension),
        120,
        "feasibility",
        false
    )
end

function all_ones_objective(prlp::PRLP)::ObjectiveFunction
    return ObjectiveFunction(
        ones(SCIP.SCIP_Real, prlp.dimension),
        10,
        "all_ones",
        false
    )
end

"""
    pstar_objective(
    p_star::Point
)::ObjectiveFunction

Construct an objective function from a point
"""

function pstar_objective(
    p_star::Point
)::ObjectiveFunction
    return ObjectiveFunction(
        as_dense_vector(p_star),
        10,
        "pstar",
        false
    )
end

function pstar_feasibility_objective(
    prlp::PRLP
)::ObjectiveFunction
    return ObjectiveFunction(
        zeros(SCIP.SCIP_Real, prlp.dimension),
        120,
        "pstar_feasibility",
        true
    )
end

"""
    point_objective(
    point::Point,
    index::Int64
)::ObjectiveFunction

Construct an objective function from a point
"""
function point_objective(
    point::Point,
    index::Int64
)::ObjectiveFunction
    return ObjectiveFunction(
        as_dense_vector(point),
        10,
        "point$(index)",
        true
    )
end

"""
    ray_objective(
    ray::Ray,
    index::Int64
)::ObjectiveFunction

Construct an objective function from a ray
"""
function ray_objective(
    ray::Ray,
    index::Int64
)::ObjectiveFunction
    return ObjectiveFunction(
        as_dense_vector(ray),
        10,
        "ray$(index)",
        true
    )
end

"""
    ObjectivePool

A datastructure to manage objective directions for PRLP
"""
mutable struct ObjectivePool
    prlp::PRLP
    point_ray_collection::PointRayCollection
    p_star_index::Int64
    remaining_indices::Vector{Int64}
end

function ObjectivePool(
    prlp::PRLP,
    point_ray_collection::PointRayCollection
)::ObjectivePool
    n_points = length(get_points(point_ray_collection))
    n_rays = length(get_rays(point_ray_collection))
    rays = get_rays(point_ray_collection)
    point_indices = collect(1:n_points)
    ray_indices = collect(1:n_rays)
    sort!(ray_indices; by = i -> abs(rays[i].objective_coefficient))
    ray_indices .+= n_points
    return ObjectivePool(prlp, point_ray_collection, -1, [point_indices; ray_indices])
end

# Utility methods for ObjectivePool
"""
    PRLPgetRowActivity(prlp::PRLP)

Get the row activity from the last PRLP solve
"""
function PRLPgetRowActivity(prlp::PRLP)
    row_activity = zeros(SCIP.SCIP_Real, length(prlp.points) + length(prlp.rays))
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(
        prlp.lpi, C_NULL, C_NULL, C_NULL, pointer(row_activity), C_NULL
    )
    return row_activity
end

function Base.iterate(op::ObjectivePool, state = 1)
    if state == 1
        return feasibility_objective(op.prlp), 2
    end
    if state == 2
        return all_ones_objective(op.prlp), 3
    end
    if state == 3
        points = get_points(op.point_ray_collection)
        op.p_star_index = argmin(1:length(points)) do i
            return get_objective_value(points[i])
        end
        return pstar_objective(points[op.p_star_index]), 4
    end
    if state == 4
        PRLPtighten(op.prlp, op.p_star_index)
        return pstar_feasibility_objective(op.prlp), 5
    end
    if state >= 5
        clean_indices(op)
        if (isempty(op.remaining_indices))
            return nothing
        end
        next_objective_index = popfirst!(op.remaining_indices)
        if next_objective_index <= length(op.prlp.points)
            point = op.prlp.points[next_objective_index]
            return (point_objective(point, next_objective_index), state + 1)
        else
            next_objective_index = next_objective_index - length(op.prlp.points)
            ray = op.prlp.rays[next_objective_index]
            return (ray_objective(ray, next_objective_index), state + 1)
        end
    end
end

function clean_indices(op::ObjectivePool)
    row_activity = PRLPgetRowActivity(op.prlp)
    filter!(op.remaining_indices) do i
        if i <= length(op.prlp.points)
            return !is_EQ(row_activity[i], 1.0)
        else
            return !is_zero(row_activity[i])
        end
    end
end
