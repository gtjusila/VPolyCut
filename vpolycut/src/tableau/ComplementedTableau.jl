using SCIP
using Lazy

# The complemented tableau is a tableau of a problem in SCIP general form 
# with the following additionals rules
# - If a column variable x is at upper bound then replace the variable with x' := -x in the tableau
# - If a row/slack variable y is at an lower bound (i.e. maximal row activity) then replace the variable with y' := -y in the tableau
# WARNING! The ConstraintMatrix does not get complemented
struct ComplementedTableau
    complemented_tableau::Tableau
    complemented_columns::Vector{Int}
end

function ComplementedTableau(tableau::Tableau)
    noriginalcols = get_noriginalcols(tableau)
    noriginalrows = get_noriginalrows(tableau)
    complemented_columns = Vector{Int}()
    @info "Original Tableau" tableau.tableau_matrix
    for i in 1:get_nvars(tableau)
        var = get_var_from_column(tableau, i)
        println(i, " Status ", get_basis_status(var))
        if is_at_upper_bound(var)
            println("Complementing $i")
            complement_column(tableau, var)
            push!(complemented_columns, i)
            println(var)
        end
    end

    return ComplementedTableau(tableau, complemented_columns)
end

function copy_and_complement(tableau::Tableau)
    new_tableau = deepcopy(tableau)
    return ComplementedTableau(new_tableau)
end

function complement_column(tableau::Tableau, var::Variable)
    set_basis_status!(var, SCIP.SCIP_BASESTAT_LOWER)
    lb = get_lb(var)
    set_lb!(var, -get_ub(var))
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
    for i in complemented_tableau.complemented_columns
        seperating_sol[i] = -seperating_sol[i]
    end
    return seperating_sol
end