using VPolyCut
using Printf

@kwdef mutable struct LPTableau <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
end

function SCIP.exec_lp(sepa::LPTableau)
    tableau = VPolyCut.construct_tableau(scip)

    print_lp_tableau(tableau)
end

function print_matrix(matrix)
    for row in eachrow(matrix)
        formatted_row = [@sprintf("%6.2f", x) for x in row]
        println(join(formatted_row, " "))
    end
end

function print_lp_tableau(lp_tableau::VPolyCut.Tableau)
    # Print the tableau matrix
    println("Tableau Matrix:")
    print_matrix(lp_tableau)

    # Print the infromation for each column
    println("================")
    println("Columns:")
    println("================")

    for i in 1:size(lp_tableau)[2]
        col = VPolyCut.get_var_from_column(lp_tableau, i)
        if isa(col, VPolyCut.LPColumn)
            println("Basis Status: ", VPolyCut.get_basis_status(col))
            println("Upper Bound: ", VPolyCut.get_ub(col))
            println("Lower Bound: ", VPolyCut.get_lb(col))
            println("Solution: ", VPolyCut.get_sol(col))
        elseif isa(col, VPolyCut.LPRow)
            println("Basis Status: ", VPolyCut.get_basis_status(col))
            println("RHS: ", VPolyCut.get_ub(col))
            println("LHS: ", VPolyCut.get_lb(col))
            println("Slack: ", VPolyCut.get_sol(col))
        end
    end

    println("================")
    println("Rows")
    println("================")

    for i in 1:VPolyCut.get_nbasis(lp_tableau)
        row = VPolyCut.get_var_from_row(lp_tableau, i)
        if isa(row, VPolyCut.LPRow)
            println("Basis Status: ", VPolyCut.get_basis_status(row))
            println("RHS: ", VPolyCut.get_ub(row))
            println("LHS: ", VPolyCut.get_lb(row))
            println("Slack: ", VPolyCut.get_sol(row))
        elseif isa(row, VPolyCut.LPColumn)
            println("Basis Status: ", VPolyCut.get_basis_status(row))
            println("Upper Bound: ", VPolyCut.get_ub(row))
            println("Lower Bound: ", VPolyCut.get_lb(row))
            println("Solution: ", VPolyCut.get_sol(row))
        end
    end

    corner = VPolyCut.construct_corner_polyhedron(lp_tableau)
    intersection_points = corner.lp_sol
    parallel_ray = corner.lp_rays
    println("====================")
    println("Intersection Points")
    println(intersection_points)
    println("Parallel Ray")
    println(parallel_ray)
    println("====================")
    return SCIP.SCIP_DIDNOTFIND
end
