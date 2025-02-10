using SCIP
function get_cut_from_separating_solution(
    tableau::Tableau, nonbasicspace::NonBasicSpace,
    separating_solution::Vector{SCIP.SCIP_Real}
)::Cut
    lp_solution = get_solution_vector(tableau)

    separating_solution = revert_point_to_original_space(nonbasicspace, separating_solution)
    b = dot(separating_solution, lp_solution) + 1

    cut_vector, b = convert_standard_inequality_to_general(
        tableau, separating_solution, b
    )

    # We normalize the cut to the form ax <= b
    return Cut(-cut_vector, -b)
end