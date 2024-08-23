using SCIP
using Random
using LinearAlgebra
using JuMP

# Constants
DEPTH_LIMIT = 100 
EPSILON = 1e-6

# Debug Flag
PRINT_OPTIMAL_TABLEAU = false


mutable struct IntersectionSeparator <: SCIP.AbstractSeparator
    called::Int64
    scipd::SCIP.SCIPData 
    row_num::Int64
    col_num::Int64
    lp_rows::Vector{Ptr{SCIP.SCIP_Row}}
    lp_cols::Vector{Ptr{SCIP.SCIP_Col}}
    IntersectionSeparator(inner) = new(0,inner,0,0,[],[])
end

function get_corner_polyhedron(sepa::IntersectionSeparator)
    """
    Get Corner Polyhedron Info From LP Solver

    For an LP min c'x s.t. Ax ≤ b, x ≥ 0 and an basic feasible solution x' and A' the extended matrix (inclusive slack).
    Let B denote the set of columns corresponding to the basic variables (inclusive slack) and N the rest of the columns
    Then rearrange A' = [B N] and x' = [x_B x_N] rearranged accordingly. The rows of the optimal tableau is hence given by
        x_i = B^{-1}b_i - B^{-1}Nx_N
    where i∈B or equivalently
        B^{-1}b = x_i + B^{-1}Nx_N = B^{-1}Ax'
    This equations hold true for any element of the LP feasible set and x' is given by setting all the components of x_N to 0
    """ 

    update_lp_data(sepa)

    col_num =  sepa.col_num
    row_num =  sepa.row_num
    ray_collection_ = []
    
    # Get Basic Indices
    basis_indices = zeros(Cint,row_num)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(sepa.scipd, pointer(basis_indices))

    #Collect Tableau Rows
    tableau = Dict()
    for (k, index) in enumerate(basis_indices)
        if index >= 0
            # Use Julia Convention of Array Indexing Starting From 1 instead of C
            tableau[index + 1] = zeros(row_num + col_num)
            SCIP.@SCIP_CALL SCIP.LibSCIP.SCIPgetLPBInvRow(sepa.scipd, k-1, pointer(tableau[index + 1],col_num+1), C_NULL, C_NULL)
            SCIP.@SCIP_CALL SCIP.LibSCIP.SCIPgetLPBInvARow(sepa.scipd, k-1, pointer(tableau[index + 1],col_num+1), pointer(tableau[index + 1]),C_NULL, C_NULL)
        end
    end

    # Get LP Columns 
    cols_ = sepa.lp_cols

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
        for j in keys(tableau)
            ray_[j] = factor_*tableau[j][i]
        end
        push!(ray_collection_, ray_)
    end

    # Generate Rays from slack variable
    rows_ = sepa.lp_rows 
    
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
        for j in keys(tableau)
            ray_[j] = factor*tableau[j][col_num+i]
            ray_non_zero |= (tableau[j][col_num+i] != 0)
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


function getSplitVariable(scipd::SCIP.SCIPData)
    choices = []
    split_variable = -1
    split_variable_pointer = nothing

    n =  SCIP.LibSCIP.SCIPgetNLPCols(scipd)
    cols = SCIP.LibSCIP.SCIPgetLPCols(scipd)
    cols = unsafe_wrap(Vector{Ptr{SCIP.LibSCIP.SCIP_COL}},cols,n)
    # Determine Splitting Variable
    for i=1:n 
        # Loop through each variable
        var = SCIP.LibSCIP.SCIPcolGetVar(cols[i])
        sol = SCIP.LibSCIP.SCIPvarGetLPSol(var)
        if SCIP.SCIPvarIsIntegral(var)==1 && (sol - floor(sol) > EPSILON)
            #Only Consider The Split if var is integral and the current solution is non integral
            push!(choices, tuple(var, i, abs(sol - (floor(sol)+0.5))))
        end
    end

    if (length(choices) == 0)
        return [] 
    end
    sort!(choices, by = x->x[end])
    split_variable = choices[1][2]
    split_variable_pointer = choices[1][1]

    return [split_variable]
end
function printCut(coef,b)
    str = ""
    for (index,i) in enumerate(coef)
        str = str*" $i"*"x$index"
    end
    str = str*" ≥ $b"
    return str 
end
function computeIntersectionPoints(scipd::SCIP.SCIPData,current_solution, split_index, ray_collection)
    intersection_points = []
    pararrel_ray = []
    for ray in ray_collection
        low = floor(current_solution[split_index])
        up = ceil(current_solution[split_index])  
        epsilon = 0 
        if ray[split_index] < 0
            epsilon = (low - current_solution[split_index])/ray[split_index]
        elseif ray[split_index] > 0 
            epsilon = (up - current_solution[split_index])/ray[split_index]
        else
            push!(pararrel_ray,ray)
        end
        push!(intersection_points, current_solution + epsilon*ray)
    end
    return intersection_points
end

function update_lp_data(sepa::IntersectionSeparator)
    """Get LP data that are important for the separator"""    
    
    row_num = Ref{Cint}(0)
    lp_rows = Ref{Ptr{Ptr{SCIP.SCIP_Row}}}(C_NULL)
    col_num = Ref{Cint}(0)
    lp_cols = Ref{Ptr{Ptr{SCIP.SCIP_Col}}}(C_NULL)
    
    SCIP.@SCIP_CALL SCIP.SCIPgetLPRowsData(sepa.scipd, lp_rows, row_num)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPColsData(sepa.scipd, lp_cols, col_num)

    sepa.row_num = row_num[]
    sepa.col_num = col_num[]
    sepa.lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}},lp_rows[],sepa.row_num)
    sepa.lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}},lp_cols[],sepa.col_num)

end

function SCIP.exec_lp(sepa::IntersectionSeparator)
    
    #Only limit 1 call 
    if sepa.called >=1
        return SCIP.SCIP_DIDNOTRUN
    end
    
    # Precondition
    if SCIP.SCIPisLPSolBasic(sepa.scipd) == false
        @warn "Solution is non basic"
        return SCIP.SCIP_DIDNOTRUN
    end
    if SCIP.SCIPgetLPSolstat(sepa.scipd) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        @warn "LP is not optimal"
        return SCIP.SCIP_DIDNOTRUN
    end
    if !(SCIP.SCIPgetStage(sepa.scipd) == SCIP.SCIP_STAGE_SOLVING)
        @warn "SCIP not solving"
        return SCIP.SCIP_DIDNOTRUN
    end

    # Write LP of current node file
    name = string(sepa.called)
    SCIP.@SCIP_CALL SCIP.SCIPwriteLP(sepa.scipd, "./temp/"*name*".lp")
    update_lp_data(sepa) 
    # Get Current LP solution 
    lp_sol, rays_collection = get_corner_polyhedron(sepa)

    
    dim = length(lp_sol)
    if length(rays_collection) != dim
        @warn "Not enough ray"
        return SCIP.SCIP_DIDNOTRUN 
    end
    
    split_var = getSplitVariable(sepa.scipd)

    if length(split_var) == 0
        @warn "No Splitting Variable"
        return SCIP.SCIP_DIDNOTRUN
    end

    split_var = split_var[1]
    sepa.called += 1

    intersection_points = computeIntersectionPoints(sepa.scipd, lp_sol, split_var, rays_collection)
    lpi = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(lpi, C_NULL, "Seperating LP", SCIP.SCIP_OBJSEN_MINIMIZE)
    # Add Variables
    SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(lpi[],dim,ones(dim),-SCIP.SCIPinfinity(sepa.scipd)*ones(dim),SCIP.SCIPinfinity(sepa.scipd)*ones(dim),C_NULL, 0,C_NULL, C_NULL, C_NULL )
  
    for point in intersection_points
        lhs = [1.0]
        rhs = [1.0]
        point = point - lp_sol
        buffer = zeros(dim)
        ind = zeros(Cint,dim)
        nnonz = 0
        j = 1
        for i in 1:dim 
            if abs(point[i]) > 10e-6
                nnonz +=1
                buffer[j] = point[i]
                ind[j] = i-1
                j+=1
            end
        end
        beg = [0]
        SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(lpi[],1,pointer(lhs),pointer(rhs),C_NULL,nnonz, pointer(beg),pointer(ind), pointer(buffer))
    end
    SCIP.@SCIP_CALL SCIP.SCIPlpiWriteLP(lpi[],"./temp/Cut$name.lp")

    # Solve The LP
    if SCIP.SCIPlpiHasPrimalSolve() == 0
        @warn "LP Solver Does Not Support Primal Solve"
        return SCIP.SCIP_DIDNOTRUN
    end
    SCIP.@SCIP_CALL SCIP.SCIPlpiSolvePrimal(lpi[])
    if SCIP.SCIPlpiIsPrimalFeasible(lpi[]) == 0
        @warn "LP Solver Failed TO Solve Cut Generation LP"
        return SCIP.SCIP_DIDNOTRUN
    end
    obj = Ref{SCIP.SCIP_Real}(0.0)
    sol = zeros(SCIP.SCIP_Real, dim)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(lpi[], obj, pointer(sol),C_NULL, C_NULL, C_NULL)
    b = dot(sol,lp_sol) + 1
    nnonz = 0
    for s in sol
        if abs(s) >= EPSILON
            nnonz +=1
        end
    end
    row = Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    cols = SCIP.LibSCIP.SCIPgetLPCols(sepa.scipd)
    cols = unsafe_wrap(Vector{Ptr{SCIP.LibSCIP.SCIP_COL}},cols,dim)
    cols__ = Array{Ptr{SCIP.SCIP_ROW}}(undef,nnonz)
    vals__ = Array{SCIP.SCIP_Real}(undef,nnonz)
    infeasible = Ref{SCIP.SCIP_Bool}()
    cnt = 1 
    for (index,s) in enumerate(sol)
        if abs(s) >= EPSILON
            cols__[cnt] = cols[index]
            vals__[cnt] = s
            cnt += 1
        end
    end
    SCIP.@SCIP_CALL SCIP.SCIPcreateRowSepa(sepa.scipd, row, sepa.scipd.sepas[sepa], "",nnonz,pointer(cols__),pointer(vals__),b,SCIP.SCIPinfinity(sepa.scipd), true,false, false)
    SCIP.@SCIP_CALL SCIP.SCIPaddRow(sepa.scipd,row[],true, infeasible)

    return SCIP.SCIP_SEPARATED
end
