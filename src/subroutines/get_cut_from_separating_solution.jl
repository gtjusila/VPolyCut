using SCIP

function get_cut_from_separating_solution(
    lp_solution::Point, nonbasicspace::NonBasicSpace,
    separating_solution::Vector{SCIP.SCIP_Real}, constraint_matrix::ConstraintMatrix
)::Cut
    separating_solution = revert_cut_vector_to_original_space(
        nonbasicspace, separating_solution
    )
    b = dot(separating_solution, lp_solution) + 1

    cut_vector, b = convert_standard_inequality_to_general(
        constraint_matrix, separating_solution, b
    )

    # We normalize the cut to the form ax <= b
    return Cut(-cut_vector, -b)
end

function convert_standard_inequality_to_general(
    constraint_matrix::ConstraintMatrix, standard_row::Vector{SCIP.SCIP_Real}, b
)
    nvars = length(standard_row)
    noriginalcols = size(constraint_matrix)[2]
    general_row = zeros(noriginalcols)

    for i in 1:nvars
        if is_zero(standard_row[i])
            continue
        end
        if i <= noriginalcols
            general_row[i] = standard_row[i]
        else
            row_index = i - noriginalcols
            for j in 1:noriginalcols
                general_row[j] -=
                    standard_row[i] * get_entry(constraint_matrix, row_index, j)
            end
            b += standard_row[i] * get_constant(constraint_matrix, row_index)
        end
    end

    return general_row, b
end