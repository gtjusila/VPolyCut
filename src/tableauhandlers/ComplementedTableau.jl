using SCIP
using Lazy

# The complemented tableau is a tableau of a problem in SCIP general form 
# with the following additionals rules
# - If a variable/slack x is at upper bound then replace the variable with x' := -x in the tableau
# WARNING! The ConstraintMatrix does not get complemented
struct ComplementedTableau
    complemented_tableau::Tableau
    complemented_columns::Vector{Int}
end

get_complemented_columns(ctableau::ComplementedTableau) = ctableau.complemented_columns

# Forwarding Functions
function get_var_from_column(ctableau::ComplementedTableau, i::Int)
    get_var_from_column(ctableau.complemented_tableau, i)
end

function get_column_from_var(ctableau::ComplementedTableau, var::Variable)
    get_column_from_var(ctableau.complemented_tableau, var)
end

@forward ComplementedTableau.complemented_tableau get_nvars
@forward ComplementedTableau.complemented_tableau create_projection_to_nonbasic_space
@forward ComplementedTableau.complemented_tableau create_trivial_projection
@forward ComplementedTableau.complemented_tableau get_solution_vector
@forward ComplementedTableau.complemented_tableau get_objective_direction
@forward ComplementedTableau.complemented_tableau get_problem_variables_pointers

function convert_standard_inequality_to_general(
    scip::SCIP.SCIPData,
    ctableau::ComplementedTableau,
    standard_row::Vector{SCIP.SCIP_Real},
    b::SCIP.SCIP_Real
)
    convert_standard_inequality_to_general(
        scip, ctableau.complemented_tableau, standard_row, b
    )
end

"""
Create a ComplementedTableau structure that wraps around an existing tableau
and modify the existing tableau accordingly
"""
function ComplementedTableau(tableau::Tableau)
    complemented_columns = Vector{Int}()

    for i in 1:get_nvars(tableau)
        var = get_var_from_column(tableau, i)
        if is_at_upper_bound(var)
            complement_column(tableau, var)
            push!(complemented_columns, i)
        end
    end

    return ComplementedTableau(tableau, complemented_columns)
end

function copy_and_complement(tableau::Tableau)
    new_tableau = deepcopy(tableau)
    return ComplementedTableau(new_tableau)
end

"""
complement the variable along with its corresponding column in the tableau
"""
function complement_column(tableau::Tableau, var::Variable)
    if get_basis_status(var) == SCIP.SCIP_BASESTAT_LOWER
        set_basis_status!(var, SCIP.SCIP_BASESTAT_UPPER)
    elseif get_basis_status(var) == SCIP.SCIP_BASESTAT_UPPER
        set_basis_status!(var, SCIP.SCIP_BASESTAT_LOWER)
    end
    lb = get_lb(var)
    set_lb!(var, -get_ub(var))
    set_obj!(var, -get_obj(var))
    set_ub!(var, -lb)
    set_sol!(var, -get_sol(var))
    col_idx = get_column_from_var(tableau, var)
    tableau[:, col_idx] = -tableau[:, col_idx]
end

function get_tableau(complemented_tableau::ComplementedTableau)
    return complemented_tableau.complemented_tableau
end

function get_uncomplemented_vector(
    seperating_sol::Vector{SCIP.SCIP_Real}, complemented_tableau::ComplementedTableau
)::Vector{SCIP.SCIP_Real}
    # Force a copy
    seperating_sol = copy(seperating_sol)
    for i in complemented_tableau.complemented_columns
        seperating_sol[i] = -seperating_sol[i]
    end
    return seperating_sol
end