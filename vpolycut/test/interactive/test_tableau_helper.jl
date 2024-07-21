import VPolyCut
using Printf

@kwdef mutable struct LPTableau <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
end

function SCIP.exec_lp(sepa::LPTableau)

    tableau = VPolyCut.get_tableau_from_scip(scip)

    print_lp_tableau(tableau)

end

function print_matrix(matrix)
    for row in eachrow(matrix)
        formatted_row = [@sprintf("%6.2f", x) for x in row]
        println(join(formatted_row, " "))
    end
end

function print_lp_tableau(lp_tableau::VPolyCut.LPTableau)
    # Print the tableau matrix
    println("Tableau Matrix:")
    print_matrix(lp_tableau)

    # Print the infromation for each column
    println("================")
    println("Columns:")
    println("================")

    for i = 1:size(lp_tableau)[2]
        col = VPolyCut.get_column_object(lp_tableau, i)
        if isa(col, VPolyCut.LPColumn)
            println("Basis Status: ", VPolyCut.get_basis_status(col))
            println("Upper Bound: ", VPolyCut.get_col_ub(col))
            println("Lower Bound: ", VPolyCut.get_col_lb(col))
            println("Solution: ", VPolyCut.get_col_sol(col))
        elseif isa(col, VPolyCut.LPRow)
            println("Basis Status: ", VPolyCut.get_basis_status(col))
            println("RHS: ", VPolyCut.get_row_rhs(col))
            println("LHS: ", VPolyCut.get_row_lhs(col))
            println("Slack: ", VPolyCut.get_row_slack(col))
        end
    end

    println("================")
    println("Rows")
    println("================")

    for i = 1:VPolyCut.get_num_lp_rows(lp_tableau)
        row = VPolyCut.get_row_object(lp_tableau, i)
        if isa(row, VPolyCut.LPRow)
            println("Basis Status: ", VPolyCut.get_basis_status(row))
            println("RHS: ", VPolyCut.get_row_rhs(row))
            println("LHS: ", VPolyCut.get_row_lhs(row))
            println("Slack: ", VPolyCut.get_row_slack(row))
        elseif isa(row, VPolyCut.LPColumn)
            println("Basis Status: ", VPolyCut.get_basis_status(row))
            println("Upper Bound: ", VPolyCut.get_col_ub(row))
            println("Lower Bound: ", VPolyCut.get_col_lb(row))
            println("Solution: ", VPolyCut.get_col_sol(row))
        end
    end

    return SCIP.SCIP_DIDNOTFIND
end