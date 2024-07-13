import VPolyCut
using Printf

@kwdef mutable struct LPTableau <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
end

function SCIP.exec_lp(sepa::LPTableau)
    SCIP.@SCIP_CALL SCIP.SCIPwriteLP(scip, "/Users/gtjusila/Documents/Project/VPolyCut/VPolyCut/test.lp")

    tableau = VPolyCut.get_scip_tableau(sepa.scipd)

    println("================")
    println("LP Tableau")
    println("================")
    print_matrix(tableau)

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