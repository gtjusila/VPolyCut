
"""
    NonBasicSpace

A NonBasicSpace object represent the non-basic space associated with an LP tableau.
It stores information about the non-basic variables, the origin point of the non basic space and which columns are complemented
"""

struct NonBasicSpace
    nonbasic_indices::Vector{Int64}
    complemented_columns::Set{Int64}
    origin_point::Point
end

function dimension(nbspace::NonBasicSpace)
    return length(nbspace.nonbasic_indices)
end

"""
    NonBasicSpace(tableau::Tableau)

Create a NonBasicSpace object from a given tableau
"""
function NonBasicSpace(tableau::Tableau)
    nonbasic_indices = Vector{Int64}()
    complemented_columns = Set{Int64}()
    mask = Vector{Int64}()
    n_tableau_vars = get_nvars(tableau)
    origin_point = get_solution_vector(tableau)
    for i in 1:n_tableau_vars
        var = get_var_from_column(tableau, i)
        if !is_basic(var)
            # Mark that the var is non-basic
            push!(nonbasic_indices, i)
            # Mark non-basic indices that must be complemented
            if is_at_upper_bound(var)
                push!(mask, length(nonbasic_indices))
            end
        end
        if is_at_upper_bound(var)
            # Mark that the var is complemented
            push!(complemented_columns, i)
        end
    end
    return NonBasicSpace(nonbasic_indices, complemented_columns, origin_point)
end

"""
    NonBasicSpace(scip::SCIP.SCIPData)

Create a NonBasicSpace object from a given SCIP object
"""
function NonBasicSpace(scip::SCIP.SCIPData)::NonBasicSpace
    basis_status = get_basis_status(scip)
    nonbasic_indices = findall(
        x -> (x == SCIP.SCIP_BASESTAT_LOWER || x == SCIP.SCIP_BASESTAT_UPPER), basis_status
    )
    complemented_columns = findall(x -> x == SCIP.SCIP_BASESTAT_UPPER, basis_status)
    origin_point = get_solution_vector(scip)
    return NonBasicSpace(nonbasic_indices, Set(complemented_columns), origin_point)
end

"""
    project_point_to_nonbasic_space(nbspace::NonBasicSpace, point::Point)

Project a point to the non-basic space
"""
function project_point_to_nonbasic_space(
    nbspace::NonBasicSpace,
    point::Point
)::Point
    # To project a point we need to remove the basic variables and complement the necessary non-basic columns
    point = point - nbspace.origin_point
    new_point = Point(dimension(nbspace))
    for i in 1:dimension(nbspace)
        idx = nbspace.nonbasic_indices[i]
        new_point[i] = ((idx in nbspace.complemented_columns) ? -1 : 1) * point[idx]
    end
    return new_point
end

"""
    revert_cut_vector_to_nonbasic_space(nbspace::NonBasicSpace, cut_vector::Vector{SCIP.SCIP_Real})

Revert a cut vector to the original space (i.e. add the basic variables and undo the complementation) 
"""
function revert_cut_vector_to_original_space(
    nbspace::NonBasicSpace,
    cut_vector::Vector{SCIP.SCIP_Real}
)::Vector{SCIP.SCIP_Real}
    new_cut_vector = zeros(SCIP.SCIP_Real, length(nbspace.origin_point))

    # Add the basic variables
    for (idx, val) in enumerate(cut_vector)
        new_cut_vector[nbspace.nonbasic_indices[idx]] = val
    end

    # Undo the complementation
    for idx in nbspace.complemented_columns
        new_cut_vector[idx] = -new_cut_vector[idx]
    end

    return new_cut_vector
end