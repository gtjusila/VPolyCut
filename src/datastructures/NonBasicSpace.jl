#
# A NonBasicSpace object represent the non-basic space associated with an LP tableau.
# It stores information about the non-basic variables, the origin point of the non basic space and which columns are complemented
# It also supports conversion of points and rays between the original space and the non-basic space
#
struct NonBasicSpace
    nonbasic_indices::Vector{Int64}
    complement_nonbasic_mask::Vector{Int64}
    complemented_columns::Set{Int64}
    origin_point::Point
end

function dimension(nbspace::NonBasicSpace)
    return length(nbspace.nonbasic_indices)
end

"""
    NonBasicSpace(scip::SCIP.SCIPData)

Create a NonBasicSpace object from a given SCIPData object.

We need to collect:
1. The LP solution which is (x,-y) where x is the solution to the problem variables and y is the row activity.
2. The complemented columns and rows which are the indices of columns which are at their upper bound and rows which are at their lower bound.
3. Finally, we should only store values in the LP solution if they are non-zero
"""
function NonBasicSpace(scip::SCIP.SCIPData)
    origin_point = get_solution_vector(scip)
    # get_basis_status directly reverse the basis status for the rows
    basis_status = get_basis_status(scip)
    nonbasic_indices = findall(
        x -> (x == SCIP.SCIP_BASESTAT_LOWER || x == SCIP.SCIP_BASESTAT_UPPER),
        basis_status
    )
    complemented_columns = Set{Int64}(
        findall(x -> x == SCIP.SCIP_BASESTAT_UPPER, basis_status)
    )
    complement_nonbasic_mask = findall(
        i -> basis_status[i] == SCIP.SCIP_BASESTAT_UPPER, nonbasic_indices
    )
    return NonBasicSpace(
        nonbasic_indices, complement_nonbasic_mask, complemented_columns, origin_point
    )
end

function project_point_to_nonbasic_space(
    nbspace::NonBasicSpace,
    point::Point
)::Point
    # To project a point we need to remove the basic variables and complement the necessary non-basic columns
    point = substract(point, nbspace.origin_point)
    new_point = Point(dimension(nbspace), get_objective_value(point))
    for i in 1:dimension(nbspace)
        idx = nbspace.nonbasic_indices[i]
        new_point[i] = ((idx in nbspace.complemented_columns) ? -1 : 1) * point[idx]
    end
    return new_point
end

function project_ray_to_nonbasic_space(
    nbspace::NonBasicSpace,
    ray::Ray
)::Ray
    # To project a ray we need to remove the basic variables and complement the necessary non-basic columns
    new_ray = Ray(ray.coefficients[nbspace.nonbasic_indices], get_generating_variable(ray))
    new_ray[nbspace.complement_nonbasic_mask] .= -new_ray[nbspace.complement_nonbasic_mask]
    return new_ray
end

function complement_ray!(nbspace::NonBasicSpace, ray::Ray)
    ray.coefficients[nbspace.complement_nonbasic_mask] .=
        -ray.coefficients[nbspace.complement_nonbasic_mask]
end
function revert_point_to_original_space(
    nbspace::NonBasicSpace,
    point::Point
)
    # To revert a point we need to add the basic variables and complement the necessary non-basic columns
    new_point = zeros(SCIP.SCIP_Real, length(nbspace.origin_point))

    for (idx, val) in enumerate(point)
        new_point[nbspace.nonbasic_indices[idx]] = val
    end

    for idx in nbspace.complemented_columns
        new_point[idx] = -new_point[idx]
    end

    return new_point
end

function revert_cut_vector_to_original_space(
    nbspace::NonBasicSpace,
    cut_vector::Vector{SCIP.SCIP_Real}
)
    # To revert a cut vector we need to add the basic variables and complement the necessary non-basic columns
    new_cut_vector = zeros(SCIP.SCIP_Real, length(nbspace.origin_point))

    for (idx, val) in enumerate(cut_vector)
        new_cut_vector[nbspace.nonbasic_indices[idx]] = val
    end

    for idx in nbspace.complemented_columns
        new_cut_vector[idx] = -new_cut_vector[idx]
    end

    return new_cut_vector
end