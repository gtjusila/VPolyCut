using SCIP

"""
    get_cut_from_separating_solution(
        separating_solution::Vector{SCIP.SCIP_Real}, constraint_matrix::ConstraintMatrix, nonbasicspace::NonBasicSpace
    )::Cut
Given a separating solution in the nonbasic space, return the corresponding cut in the original space
"""
function get_cut_from_separating_solution(
    separating_solution::Vector{SCIP.SCIP_Real},
    nonbasicspace::NonBasicSpace,
    beta::SCIP.SCIP_Real)::Cut
    origin_point = nonbasicspace.origin_point
    constraint_matrix = nonbasicspace.constraint_matrix

    separating_solution = revert_cut_vector_to_original_space(
        nonbasicspace, separating_solution
    )
    b = dot(separating_solution, origin_point) + beta
    n_rows, n_cols = size(constraint_matrix)

    for i in 1:n_rows
        if !is_zero(separating_solution[i + n_cols])
            # Substitute slack with original variables

            # We need to find the non-zero entries in the constraint row corresponding to the slack
            indices, values = findnz(constraint_matrix.data[i, :])

            # If the inequality is a^Tx + y.s <= b, and u^Tx + s + c = 0 then the inequality can be 
            # rewritten as (a-y.u)^Tx <= b + y.c
            for (j, val) in zip(indices, values)
                separating_solution[j] -= val * separating_solution[i + n_cols]
            end
            b += separating_solution[i + n_cols] * constraint_matrix.constants[i]
        end
    end

    # We can now remove the slack variable from the separating solution 
    separating_solution = separating_solution[1:n_cols]
    # We normalize the cut to the form ax <= b
    return Cut(-separating_solution, -b)
end
