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
    On fail exception
    """
    on_fail::Union{Exception,Nothing}
    """
    label
    """
    label::String
end

# Constructor methods for different objective functions

function feasibility_objective(prlp::PRLP)::ObjectiveFunction
    return ObjectiveFunction(
        zeros(SCIP.SCIP_Real, prlp.dimension),
        120,
        FailedToProvePRLPFeasibility(),
        "feasibility"
    )
end

function all_ones_objective(prlp::PRLP)::ObjectiveFunction
    return ObjectiveFunction(
        ones(SCIP.SCIP_Real, prlp.dimension),
        10,
        nothing,
        "all_ones"
    )
end

function pstar_objective(
    p_star::Point
)::ObjectiveFunction
    return ObjectiveFunction(
        as_dense_vector(p_star),
        10,
        nothing,
        "pstar"
    )
end

function pstar_feasibility_objective(
    prlp::PRLP
)::ObjectiveFunction
    return ObjectiveFunction(
        zeros(SCIP.SCIP_Real, prlp.dimension),
        120,
        PStarNotTight(),
        "pstar_feasibility"
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
        next_objective = popfirst!(op.remaining_indices)
        is_point = next_objective <= length(op.prlp.points)
        if is_point
            return (
                ObjectiveFunction(
                    as_dense_vector(op.prlp.points[next_objective]),
                    10,
                    nothing,
                    "point_$(next_objective)"
                ),
                state + 1)
        else
            return (
                ObjectiveFunction(
                    as_dense_vector(op.prlp.rays[next_objective - length(op.prlp.points)]),
                    10,
                    nothing,
                    "ray_$(next_objective - length(op.prlp.points))"
                ),
                state + 1)
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
