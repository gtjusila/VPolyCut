# 
# This file contains the code to retrieve the corner polyhedron from the current SCIP Pointer 
# along with the helper function to get Tableau Information
#
import SCIP

function get_corner_polyhedron(scip::SCIP.SCIPData)
    """
    Get Corner Polyhedron Information From LP Solver from the given scip pointer
    """
    
    # C Style Allocation
    row_num = Ref{Cint}(0)
    lp_rows = Ref{Ptr{Ptr{SCIP.SCIP_Row}}}(C_NULL)
    col_num = Ref{Cint}(0)
    lp_cols = Ref{Ptr{Ptr{SCIP.SCIP_Col}}}(C_NULL)
    
    # Get Necessary SCIP Data 
    SCIP.@SCIP_CALL SCIP.SCIPgetLPRowsData(scip, lp_rows, row_num)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPColsData(scip, lp_cols, col_num)

    # Reformat into julia style data type
    row_num = row_num[]
    col_num = col_num[]
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}},lp_rows[],row_num)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}},lp_cols[],col_num)

    # Initiate a vector to collect corner polyhedron ray
    ray_collection_ = Vector{Vector{Int64}}(undef,0) 
    
    # Determine Basic Indices
    basis_indices_ = zeros(Cint,row_num)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, pointer(basis_indices_))
    get_dense_tableau_rows(scip); 
    #Collect Tableau Rows
    tableau_ = Dict()
    for (k, index) in enumerate(basis_indices_)
        if index >= 0
            tableau_[index + 1] = zeros(row_num + col_num)
            SCIP.@SCIP_CALL SCIP.LibSCIP.SCIPgetLPBInvRow(scip, k-1, pointer(tableau_[index + 1],col_num+1), C_NULL, C_NULL)
            SCIP.@SCIP_CALL SCIP.LibSCIP.SCIPgetLPBInvARow(scip, k-1, pointer(tableau_[index + 1],col_num+1), pointer(tableau_[index + 1]),C_NULL, C_NULL)
        end
    end

    # Get LP Columns 
    cols_ = lp_cols

    #=
    Generate Rays from non basic variable
    If variable is non basic i.e. the upper or lower bound is attained then one can push the variable against the bound
    =#
    for i=1:col_num 
        # SCIP Basic Vartiables
        if (SCIP.SCIPcolGetBasisStatus(cols_[i]) == SCIP.SCIP_BASESTAT_UPPER)
            factor_ = 1.0; 
        elseif (SCIP.SCIPcolGetBasisStatus(cols_[i]) == SCIP.SCIP_BASESTAT_LOWER)
            factor_ = -1.0;
        elseif (SCIP.SCIPcolGetBasisStatus(cols_[i]) == SCIP.SCIP_BASESTAT_ZERO)
            #Safekeeping: Should Never Happen 
            return SCIP.SCIP_DIDNOTRUN
        else
            continue
        end

        ray_ = zeros(col_num)
        ray_[i] = -factor_
        for j in keys(tableau_)
            ray_[j] = factor_*tableau_[j][i]
        end
        push!(ray_collection_, ray_)
    end

    # Generate Rays from slack variable
    rows_ = lp_rows 
    
    for i=1:row_num
        #Go through each row
        if( SCIP.LibSCIP.SCIProwGetBasisStatus(rows_[i]) == SCIP.SCIP_BASESTAT_LOWER)
            factor = 1.0
        elseif (SCIP.LibSCIP.SCIProwGetBasisStatus(rows_[i]) == SCIP.SCIP_BASESTAT_UPPER)
            factor = - 1.0
        elseif (SCIP.LibSCIP.SCIProwGetBasisStatus(rows_[i]) ==  SCIP.SCIP_BASESTAT_ZERO)
            return SCIP.SCIP_DIDNOTRUN
        else
            continue
        end

        ray_ = zeros(col_num)
        ray_non_zero = false
        for j in keys(tableau_)
            ray_[j] = factor*tableau_[j][col_num+i]
            ray_non_zero |= (tableau_[j][col_num+i] != 0)
        end

        if ray_non_zero
            push!(ray_collection_, ray_)
        end
    end

    current_solution = zeros(col_num)
    for i=1:col_num
        var = SCIP.LibSCIP.SCIPcolGetVar(cols_[i])
        current_solution[i] = SCIP.LibSCIP.SCIPvarGetLPSol(var)
    end

    return current_solution,ray_collection_
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
        SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvARow(scip,i-1,C_NULL,row_,C_NULL, C_NULL)
        tableau[i] = row_
    end
    println(string(tableau))
    return
end
