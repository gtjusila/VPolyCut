# 
# This file contains the code to retrieve the corner polyhedron from the current SCIP Pointer 
# along with the helper function to get Tableau Information
#
import SCIP
include("utils.jl")

function get_corner_polyhedron(scip::SCIP.SCIPData)
    """
    Get Corner Polyhedron Information From LP Solver from the given scip pointer
    """
    @assert SCIP.SCIPgetLPSolstat(scip) == SCIP.SCIP_LPSOLSTAT_OPTIMAL
    # Get Information from SCIP
    lp_cols = get_lp_columns(scip)
    col_num = length(lp_cols)
    lp_rows,row_num = get_lp_row_information(scip) 
    basis_indices = get_lp_basis_information(scip)

    # Initiate a vector to collect corner polyhedron ray
    ray_collection = Vector{Vector{SCIP.SCIP_Real}}(undef,0) 
    
    # Get Tableau
    tableau = get_dense_tableau_rows(scip) 

    #=
    Generate Rays from non basic variable
    If variable is non basic i.e. the upper or lower bound is attained then one can push the variable against the bound
    =#
    for (j, col) in enumerate(lp_cols) 
        # SCIP Basic Vartiables
        if (SCIP.SCIPcolGetBasisStatus(col) == SCIP.SCIP_BASESTAT_UPPER)
            factor_ = 1.0; 
        elseif (SCIP.SCIPcolGetBasisStatus(col) == SCIP.SCIP_BASESTAT_LOWER)
            factor_ = -1.0;
        elseif (SCIP.SCIPcolGetBasisStatus(col) == SCIP.SCIP_BASESTAT_ZERO)
            #Safekeeping: Should Never Happen unless polyhedron is not pointed
            throw("SCIP_BASESTAT_ZERO encountered") 
        elseif (SCIP.SCIPcolGetBasisStatus(col) == SCIP.SCIP_BASESTAT_BASIC)
            # variable is in the basis 
            continue
        else
            throw("SCIP_BASESTAT Status Undefined")
        end

        # Construct ray r_j
        ray_ = zeros(col_num)
        ray_[j] = -factor_
        for (k, index) in enumerate(basis_indices)
            if index >= 0
                ray_[index+1] =  factor_*tableau[k][j]
            end
        end
        
        push!(ray_collection, ray_)
    end

    # Generate Rays from slack variable
    for (i, row) in enumerate(lp_rows) 
        #Go through each row
        if( SCIP.LibSCIP.SCIProwGetBasisStatus(row) == SCIP.SCIP_BASESTAT_LOWER)
            factor_ = 1.0
        elseif (SCIP.LibSCIP.SCIProwGetBasisStatus(row) == SCIP.SCIP_BASESTAT_UPPER)
            factor_ = - 1.0
        elseif (SCIP.LibSCIP.SCIProwGetBasisStatus(row) ==  SCIP.SCIP_BASESTAT_ZERO)
            # Same as above
            throw("SCIP_BASESTAT_ZERO encountered")
        elseif (SCIP.LibSCIP.SCIProwGetBasisStatus(row) == SCIP.SCIP_BASESTAT_BASIC)
            # Row is basic
            continue
        else
            throw("SCIP_BASESTAT Status Undefined")
        end

        # Construct ray r_j
        ray_ = zeros(col_num)
        ray_non_zero = false
        for (j, index) in enumerate(basis_indices)
            if (index >= 0) && (tableau[j][col_num+i] != 0)
                ray_[index+1] =  factor_*tableau[j][col_num+i]
                ray_non_zero =  true
            end
        end

        if ray_non_zero
            push!(ray_collection, ray_)
        end
    end

    current_solution = get_lp_solution_vector(scip; col_num = col_num, lp_cols = lp_cols) 

    return current_solution,ray_collection
end
