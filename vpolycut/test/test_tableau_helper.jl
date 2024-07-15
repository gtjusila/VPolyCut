import VPolyCut
using Printf

@kwdef mutable struct LPTableau <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
end

function SCIP.exec_lp(sepa::LPTableau)
    SCIP.@SCIP_CALL SCIP.SCIPwriteLP(scip, "/Users/gtjusila/Documents/Project/VPolyCut/VPolyCut/test.lp")

    tableau = VPolyCut.get_dense_tableau_from_scip(scip)

    print_dense_tableau_info(tableau)

    row_count = VPolyCut.get_lp_row_count(scip)
    buffer = SCIP.SCIPgetLPRows(scip)
    rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, buffer, row_count)
    #=
    println("================")
    println("Column Information")
    println("================")
    lp_cols = get_lp_columns(sepa.scipd)
    for (i, col) in enumerate(lp_cols)
        var = SCIP.SCIPcolGetVar(col)
        name = unsafe_string(SCIP.SCIPvarGetName(var))
        status = SCIP.SCIPcolGetBasisStatus(col)
        println("Column $i Variable $(name) Status $(status)")
    end

    println("================")
    println("Row Information")
    println("================")
    lp_rows = get_lp_rows(sepa.scipd)
    for (i, row) in enumerate(lp_rows)
        cons = SCIP.SCIProwGetOriginCons(row)
        name = unsafe_string(SCIP.SCIPconsGetName(cons))
        status = SCIP.SCIProwGetBasisStatus(row)
        println("Column $i Constraint $(name) Status $(status)")
    end

    println("================")
    println("LP Solution Vector")
    println("================")
    sol = get_lp_solution_vector(scip)
    sol = [round(x, sigdigits=2) for x in sol]
    println(sol)

    println("===============")
    =#
    return SCIP.SCIP_DIDNOTFIND
end

function print_matrix(matrix)
    for row in eachrow(matrix)
        formatted_row = [@sprintf("%6.2f", x) for x in row]
        println(join(formatted_row, " "))
    end
end

function get_lp_columns(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Col}}
    col_num = Ref{Cint}(0)
    lp_cols = Ref{Ptr{Ptr{SCIP.SCIP_Col}}}(C_NULL)

    SCIP.@SCIP_CALL SCIP.SCIPgetLPColsData(scip, lp_cols, col_num)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, lp_cols[], col_num[])

    return lp_cols
end

function get_lp_rows(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Row}}
    lp_rows = Ref{Ptr{Ptr{SCIP.SCIP_Row}}}(C_NULL)
    row_num = Ref{Cint}(0)
    # Get Necessary SCIP Data 
    SCIP.@SCIP_CALL SCIP.SCIPgetLPRowsData(scip, lp_rows, row_num)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, lp_rows[], row_num[])

    return lp_rows
end

function get_lp_solution_vector(scip)
    # If data is not given then first fetch data from scip
    lp_cols = get_lp_columns(scip)
    col_num = length(lp_cols)

    # Get the solution
    current_solution = zeros(col_num)
    for i = 1:col_num
        var_ = SCIP.LibSCIP.SCIPcolGetVar(lp_cols[i])
        current_solution[i] = SCIP.LibSCIP.SCIPvarGetLPSol(var_)
    end

    return current_solution
end

#
# CHATGPT CODE COMING
#
function print_dense_tableau_info(tableau::VPolyCut.DenseTableau)
    println("LP Tableau Matrix:")

    print_matrix(tableau.tableau_matrix)
    num_rows, num_columns = size(tableau.tableau_matrix)
    println("\nNumber of rows: ", num_rows)
    println("Number of columns: ", num_columns)

    println("\nColumn Information:")
    for col_idx in 1:num_columns
        col_found = false
        for col in VPolyCut.get_tableau_columns(tableau)
            if VPolyCut.get_matrix_column_index_for_tableau_column(col) == col_idx
                println("Column $col_idx belongs to an LPColumn with index $(col.index) having lower bound $(VPolyCut.get_tableau_column_lower_bound(col)) and upper bound $(VPolyCut.get_tableau_column_upper_bound(col)) and basis status $(VPolyCut.get_tableau_column_basis_status(col))")
                col_found = true
                break
            end
        end
        if !col_found
            for row in VPolyCut.get_tableau_rows(tableau)
                if VPolyCut.get_matrix_column_index_for_tableau_row(row) == col_idx
                    println("Column $col_idx belongs to an LPRow with index $(row.index) having lhs $(VPolyCut.get_tableau_row_lhs(row)) and rhs $(VPolyCut.get_tableau_row_rhs(row)) and basis status $(VPolyCut.get_tableau_row_basis_status(row))")
                    col_found = true
                    break
                end
            end
        end
        if !col_found
            println("Column $col_idx does not belong to any TableauColumn or TableauRow")
        end
    end

    println("\nRow Information:")
    for row_idx in 1:num_rows
        row_found = false
        for col in VPolyCut.get_tableau_columns(tableau)
            if VPolyCut.get_matrix_row_index_for_tableau_column(col) == row_idx
                println("Row $row_idx belongs to an LPColumn with index $(col.index) having lower bound $(VPolyCut.get_tableau_column_lower_bound(col)) and upper bound $(VPolyCut.get_tableau_column_upper_bound(col))")
                row_found = true
                break
            end
        end
        if !row_found
            for row in VPolyCut.get_tableau_rows(tableau)
                if VPolyCut.get_matrix_row_index_for_tableau_row(row) == row_idx
                    println("Row $row_idx belongs to an LPRow with index $(row.index) having lhs $(VPolyCut.get_tableau_row_lhs(row)) and rhs $(VPolyCut.get_tableau_row_rhs(row))")
                    row_found = true
                    break
                end
            end
        end
        if !row_found
            println("Row $row_idx does not belong to any TableauColumn or TableauRow")
        end
    end
end