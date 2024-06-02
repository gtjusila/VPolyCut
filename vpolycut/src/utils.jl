# Some Common Functions Needed By All Function
using Random

function get_lp_columns(scip::SCIP.SCIPData)
    
    col_num = Ref{Cint}(0)
    lp_cols = Ref{Ptr{Ptr{SCIP.SCIP_Col}}}(C_NULL)
    
    SCIP.@SCIP_CALL SCIP.SCIPgetLPColsData(scip, lp_cols, col_num)
    
    col_num = col_num[]
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}},lp_cols[],col_num)

    return lp_cols
end

function get_lp_row_information(scip::SCIP.SCIPData)
    
    row_num = Ref{Cint}(0)
    lp_rows = Ref{Ptr{Ptr{SCIP.SCIP_Row}}}(C_NULL)
    
    # Get Necessary SCIP Data 
    SCIP.@SCIP_CALL SCIP.SCIPgetLPRowsData(scip, lp_rows, row_num)

    row_num = row_num[]
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}},lp_rows[],row_num)

    return lp_rows,row_num
end

function get_lp_solution_vector(scip; col_num::Int64 = Int64(0), lp_cols::Vector{Ptr{SCIP.SCIP_Col}} = Vector{Ptr{SCIP.SCIP_Col}}(undef,0))
    # If data is not given then first fetch data from scip
    col_num = Int32(col_num)
    if col_num == 0 || length(lp_cols) == 0
        lp_cols = get_lp_columns(scip) 
        col_num = length(lp_cols)
    end
    
    # Get the solution
    current_solution = zeros(col_num)
    for i=1:col_num
        var_ = SCIP.LibSCIP.SCIPcolGetVar(lp_cols[i])
        current_solution[i] = SCIP.LibSCIP.SCIPvarGetLPSol(var_)
    end

    return current_solution
end

function get_lp_basis_information(scip)
    row_num = SCIP.SCIPgetNLPRows(scip)
    basis_indices = zeros(Cint,row_num)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, pointer(basis_indices))
    return basis_indices
end

function print_col_status(scip)
    lp_cols = get_lp_columns(scip)
    col_num = length(lp_cols)
    for (i,col) in enumerate(lp_cols)
        status = SCIP.SCIPcolGetBasisStatus(col)
        reduced_cost = SCIP.SCIPgetColRedcost(scip,col)
        println("Column "*string(i)*" status "*string(status)*" reduced cost "* string(reduced_cost))
    end
end

function print_row_status(scip)
    lp_rows, row_num = get_lp_row_information(scip)
    for (i,row) in enumerate(lp_rows)
        status = SCIP.SCIProwGetBasisStatus(row)
        println("Row "*string(i)*" status "*string(status))
    end
end

# [CLEAN] Instead Of Getting All of the tableau rows we can also just get the tableau rows corresponding
# to basic non slack variable. But for debugging and understanding purposes we get all tableau rows here
# in production this should be cleaned
function get_dense_tableau_rows(scip)

    tableau = Dict{Int64, Vector{SCIP.SCIP_Real}}()
    m = SCIP.SCIPgetNLPRows(scip) # Tableau will have M rows
    n = SCIP.SCIPgetNLPCols(scip) # LP have N Columns

    for i=1:m
        # Allocate Memory to store tableau entries
        row_ = zeros(SCIP.SCIP_Real,m+n)
        # Tableau is given by [B^{-1}A B^{-1}]
        SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvRow(scip,i-1,pointer(row_,n+1),C_NULL, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvARow(scip,i-1,pointer(row_,n+1),pointer(row_),C_NULL, C_NULL)
        tableau[i] = row_
    end

    return tableau
end

function get_lp_dual_solution(scip)
    
    lp_row, row_num = get_lp_row_information(scip)
    dual_solution = zeros(SCIP.SCIP_Real,row_num)
    
    for (i, row) in enumerate(lp_row)
        dual_solution[i] = SCIP.SCIProwGetDualsol(row)
    end
    
    return dual_solution
end

"""
The function receive as input an scip pointer and return a vector of pointers to variable currently in the LP
"""
function get_lp_variables(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Var}}
    vars = Ptr{SCIP.SCIP_Var}[]

    cols = get_lp_columns(scip)
    for col in cols
        push!(vars,SCIP.SCIPcolGetVar(col))
    end

    return vars
end
function print_lp_information(scip)
    
    println("====================")
    println("Printing LP Information\n") 

    println("LP Tableau")
    println("----------")
    println(get_dense_tableau_rows(scip))
   
    println()
    println("Is dual feasible? "*string(SCIP.SCIPisLPDualReliable(scip)))
    println("Basic Indices "*string(get_lp_basis_information(scip)))
    
    println()
    println("Row Status")
    println("----------")
    print_row_status(scip)

    println() 
    println("Col Status")
    println("----------")
    print_col_status(scip)

    println()
    println("----------")
    println("Extra information From LPI")
    println("----------")

    lpi = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPI(scip, lpi)
    primal = zeros(SCIP.SCIP_Real, 10)
    dual =  zeros(SCIP.SCIP_Real, 10)
    activity =  zeros(SCIP.SCIP_Real, 10)
    red_cost = zeros(SCIP.SCIP_Real,10)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(lpi[],C_NULL, primal, dual, activity, red_cost)
    
    println()
    println("Primal Solution")
    println("----------")
    println(primal)

    println()
    println("Dual Solution")
    println("----------")
    println(dual)
    
    println()
    println("Row Activity")
    println("----------")
    println(activity)
    
    println()
    println("Reduced Cost")
    println("----------")
    println(red_cost)
    
    name = randstring(3)
    SCIP.SCIPlpiWriteLP(lpi[],joinpath(WRITE_PATH, name*".lp"))
    println()
    println("LP written as "*name*".lp")
    println("====================")
end 

"""
Given a vector of SCIP variable returns the index of an integer variable with highest pseudobranch score
"""
function get_first_fractional_index(vars::Vector{Ptr{SCIP.SCIP_Var}})
    
    # Initialize
    split_index = -1

    # Determine Splitting Variable
    for (i, var) in enumerate(vars) 
        # Loop through each variable
        sol = SCIP.LibSCIP.SCIPvarGetLPSol(var)
        if SCIP.SCIPvarIsIntegral(var)==1 && (sol - (floor(sol)) > 0.01) && (ceil(sol) - sol > 0.01)
            #Only Consider The Split if var is integral and the current solution is non integral
            split_index = i
        end
    end

    return split_index
end
