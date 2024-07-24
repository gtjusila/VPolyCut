using SCIP
@kwdef mutable struct Projection
    original_dimension::Int = 0
    projected_dimension::Int = 0
    map_projected_to_original::Vector{Int} = []
end

function create_projection_to_nonbasic_space(tableau::Tableau)::Projection
    original_dimension = get_nvars(tableau)
    projected_to_original = []

    j = 1
    for i = 1:original_dimension
        var = get_var_from_column(tableau, i)
        if !is_basic(var)
            push!(projected_to_original, i)
            j += 1
        end
    end

    return Projection(original_dimension, j - 1, projected_to_original)
end

function project_point(projection::Projection, point::Point)
    return point[projection.map_projected_to_original]
end

function undo_projection(projection::Projection, point::Point)
    original_point = zeros(projection.original_dimension)
    for i = 1:length(projection.map_projected_to_original)
        original_point[projection.map_projected_to_original[i]] = point[i]
    end
    return original_point
end