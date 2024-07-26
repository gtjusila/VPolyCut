using SCIP
@kwdef mutable struct Projection
    original_dimension::Int = 0
    projected_dimension::Int = 0
    map_projected_to_original::Vector{Int} = []
end

"""
The nonbasic projection is the projection into the subspace spanned by the nonbasic variables
"""
function create_nonbasic_projection(tableau::Tableau)::Projection
    original_dimension = get_nvars(tableau)
    projected_to_original = []

    j = 1
    for i in 1:original_dimension
        var = get_var_from_column(tableau, i)
        if !is_basic(var)
            push!(projected_to_original, i)
            j += 1
        end
    end

    return Projection(original_dimension, j - 1, projected_to_original)
end

"""
The trivial projeciton is the projection into the subspace spanned by the original variables
"""
function create_trivial_projection(tableau::Tableau)::Projection
    original_dimension = get_nvars(tableau)
    projected_dimension = get_noriginalcols(tableau)
    map_projected_to_original = collect(1:projected_dimension)
    return Projection(original_dimension, projected_dimension, map_projected_to_original)
end
function project_point(projection::Projection, point::Point)
    return point[projection.map_projected_to_original]
end

function undo_projection(projection::Projection, point::Point)
    original_point = zeros(projection.original_dimension)
    for i in 1:length(projection.map_projected_to_original)
        original_point[projection.map_projected_to_original[i]] = point[i]
    end
    return original_point
end
