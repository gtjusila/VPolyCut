#
# A NonBasicSpace object represent the non-basic space associated with an LP tableau.
# It stores information about the non-basic variables, the origin point of the non basic space and which columns are complemented
# It also supports conversion of points and rays between the original space and the non-basic space
#
struct NonBasicSpace
    nonbasic_indices::Vector{Int64}
    complemented_columns::Set{Int64}
    origin_point::Point
end

"""
    NonBasicSpace(tableau::Tableau)

Create a NonBasicSpace object from a given tableau
"""
function NonBasicSpace(tableau::Tableau)
    nonbasic_indices = Vector{Int64}()
    complemented_columns = Set{Int64}()
    n_tableau_vars = get_nvars(tableau)
    origin_point = get_solution_vector(tableau)
    j = 1
    for i in 1:n_tableau_vars
        var = get_var_from_column(tableau, i)
        if !is_basic(var)
            # Mark that the var is non-basic
            push!(nonbasic_indices, i)
        end
        if is_at_upper_bound(var)
            # Mark that the var is complemented
            push!(complemented_columns, i)
            origin_point[j] = -origin_point[j]
        end
    end

    return NonBasicSpace(nonbasic_indices, complemented_columns, origin_point)
end

function project_point_to_nonbasic_space(
    nbspace::NonBasicSpace,
    point::Point
)
    # To project a point we need to remove the basic variables and complement the necessary non-basic columns
    point = point - nbspace.origin_point
    new_point = Point()
    for idx in nbspace.nonbasic_indices
        push!(new_point, ((idx in nbspace.complemented_columns) ? -1 : 1) * point[idx])
    end
    return new_point
end

function project_ray_to_nonbasic_space(
    nbspace::NonBasicSpace,
    ray::Ray
)
    # To project a ray we need to remove the basic variables and complement the necessary non-basic columns
    new_ray = Vector{SCIP.SCIP_Real}(undef, 0)
    for idx in nbspace.nonbasic_indices
        push!(new_ray, ((idx in nbspace.complemented_columns) ? -1 : 1) * ray[idx])
    end
    return Ray(new_ray, get_generating_variable(ray))
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