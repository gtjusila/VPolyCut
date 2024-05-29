import SCIP
using LinearAlgebra
using JuMP

CALL_LIMIT = 100
EPSILON = 1e-10
WRITE_PATH = joinpath(pwd(),"temp")

include("../src/CornerPolyhedron.jl")
include("../src/utils.jl")

DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON = false 
DEBUG_PRINT_INTERSECTION_POINTS = false

@kwdef mutable struct IntersectionSeparator <: SCIP.AbstractSeparator
    called::Int64 = 0
    debug_sol_path::String = "" 
    scipd::SCIP.SCIPData 
end

function decideSplitIndex(scipd::SCIP.SCIPData)
    
    # Initialize
    split_index = -1

    # Get Relevant SCIP Data
    n =  SCIP.LibSCIP.SCIPgetNLPCols(scipd)
    cols = SCIP.LibSCIP.SCIPgetLPCols(scipd)
    cols = unsafe_wrap(Vector{Ptr{SCIP.LibSCIP.SCIP_COL}},cols,n)
    
    # Determine Splitting Variable
    for i=1:n 
        # Loop through each variable
        var = SCIP.LibSCIP.SCIPcolGetVar(cols[i])
        sol = SCIP.LibSCIP.SCIPvarGetLPSol(var)
        if SCIP.SCIPvarIsIntegral(var)==1 && (sol - (floor(sol)) > 0.2) && (ceil(sol) - sol > 0.2)
            #Only Consider The Split if var is integral and the current solution is non integral
            split_index = i
        end
    end

    return split_index
end

function computeIntersectionPoints(current_solution, split_index, ray_collection)
    
    intersection_points = []
    pararrel_ray = []
    
    for ray in ray_collection
        low = floor(current_solution[split_index])
        up = ceil(current_solution[split_index])  
        epsilon = 0 
        if ray[split_index] < -EPSILON
            epsilon = (low - current_solution[split_index])/ray[split_index]
        elseif ray[split_index] > EPSILON
            epsilon = (up - current_solution[split_index])/ray[split_index]
        else
            push!(pararrel_ray,ray)
            continue;
        end
        push!(intersection_points, current_solution + epsilon*ray)
    end
    return (intersection_points,pararrel_ray)
end

function constructSeperatingLP(sepa, dim, lp_sol, intersection_points, parallel_ray)
    
    lpi = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(lpi, C_NULL, "Seperating LP", SCIP.SCIP_OBJSEN_MINIMIZE)

    # Add Variables
    SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(lpi[],dim,zeros(dim),-10000*ones(dim),10000*ones(dim),C_NULL, 0,C_NULL, C_NULL, C_NULL )
    for (idx,point) in enumerate(intersection_points)
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
        if nnonz >0
            SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(lpi[],1,pointer(lhs),pointer(rhs),C_NULL,nnonz, pointer(beg),pointer(ind), pointer(buffer))
        end
    end

    for ray in parallel_ray
        lhs = [0.0001]
        rhs = [SCIP.SCIPinfinity(sepa.scipd)]
        buffer = zeros(dim)
        ind = zeros(Cint,dim)
        nnonz = 0
        j = 1
        for i in 1:dim 
            if abs(ray[i]) > 10e-6
                nnonz +=1
                buffer[j] = ray[i]
                ind[j] = i-1
                j+=1
            end
        end
        beg = [0]
        if nnonz >0
            SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(lpi[],1,pointer(lhs),pointer(rhs),C_NULL,nnonz, pointer(beg),pointer(ind), pointer(buffer))
        end 
    end

    SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(lpi[],dim,ones(dim),-10000*ones(dim),10000*ones(dim),C_NULL, 0,C_NULL, C_NULL, C_NULL)
    for i=1:dim
        lhs = [-SCIP.SCIPinfinity(sepa.scipd)]
        rhs = [0.0]
        nnonz = 2
        beg = [0]
        ind = zeros(Cint,2*dim)
        buffer = zeros(2*dim)
        ind[1] = i-1
        ind[2] =  i-1+dim
        buffer[1] = -1
        buffer[2] = -1
        SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(lpi[],1,pointer(lhs),pointer(rhs),C_NULL,nnonz,pointer(beg),pointer(ind),pointer(buffer)) 
        buffer[1] = 1
        buffer[2] = -1
        SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(lpi[],1,pointer(lhs),pointer(rhs),C_NULL,nnonz,pointer(beg),pointer(ind),pointer(buffer)) 
    end

    return lpi
end

#[CLEAN] Debug Function

function check_in_corner_polyhedron(vertex, rays, point)
    model = Model(SCIP.Optimizer)
    set_attribute(model, "display/verblevel",0)
    @variable(model, l[i=1:length(rays)]>=-0.1)
    @variable(model, z)
    @constraint(model, l .>= z)
    @objective(model, Max, z)
    R = hcat(rays...)
    @constraint(model, R*l + vertex .== point)
    optimize!(model)
    if !is_solved_and_feasible(model)
        @warn "Reference Solution Not In Corner Polyhedron"
    else
        println("Solution Is In The Corner Polyhedron")
    end
end

function verify_reference_solution_in_disjunction(sepa::IntersectionSeparator,reference_solution, split_index)
    lp_cols = get_lp_columns(sepa.scipd)
    col_num = length(lp_cols)
    split_var = SCIP.SCIPcolGetVar(lp_cols[split_index])
    split_var_value = SCIP.SCIPvarGetSol(split_var,1)
    reference = SCIP.SCIPgetSolVal(sepa.scipd, reference_solution[], split_var)
    println(string(split_var_value)*" vs "*string(reference))
    
    found = false
    SCIP.@SCIP_CALL SCIP.SCIPstartProbing(sepa.scipd)
    SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(sepa.scipd)
    SCIP.@SCIP_CALL SCIP.SCIPchgVarUbProbing(sepa.scipd, split_var, floor(split_var_value))
    println("Testing For Lower Than")
    feasible = Ref{SCIP.SCIP_Bool}(0)
    SCIP.@SCIP_CALL SCIP.SCIPcheckSol(sepa.scipd,reference_solution[], 1,1,1,0,1,feasible)
    println(feasible[]==1)
    found = found || (feasible[]==1)
    SCIP.@SCIP_CALL SCIP.SCIPbacktrackProbing(sepa.scipd, 0)
    SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(sepa.scipd)
    SCIP.@SCIP_CALL SCIP.SCIPchgVarLbProbing(sepa.scipd,split_var, ceil(split_var_value))
    println("Testing For Larger Than")
    feasible = Ref{SCIP.SCIP_Bool}(0)
    SCIP.@SCIP_CALL SCIP.SCIPcheckSol(sepa.scipd,reference_solution[], 1,1,1,0,1,feasible)
    println(feasible[]==1)
    found = found || (feasible[]==1) 
    SCIP.@SCIP_CALL SCIP.SCIPendProbing(sepa.scipd)
end

function SCIP.exec_lp(sepa::IntersectionSeparator)
    
    # Check Preconditions
    @assert(SCIP.SCIPgetStage(sepa.scipd) != sepa.scipd)
    @assert(SCIP.SCIPisLPSolBasic(sepa.scipd) != 0)
    @assert(SCIP.SCIPgetLPSolstat(sepa.scipd) == SCIP.SCIP_LPSOLSTAT_OPTIMAL) 

    # Only limit 1 call 
    if sepa.called >= CALL_LIMIT
        return SCIP.SCIP_DIDNOTRUN
    end
    sepa.called +=1 
    
    # STEP 1: Get Corner Polyhedron of current LP solution
    lp_sol = nothing
    lp_rays = nothing
    lp_sol, lp_rays = get_corner_polyhedron(sepa.scipd)
    
    # [CLEAN] For Debug Write LP of current node file THIS STEP MUST BE DONE AFTER getting the corner polyhedron since otherwise LP SOLSTAT will be affected
    name = string(sepa.called)
    SCIP.@SCIP_CALL SCIP.SCIPwriteMIP(sepa.scipd, joinpath(WRITE_PATH,name*".mip"),true, false, true)

    if DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON
        println("====================")
        println("Seperator Called")
        println("LP Solution is "*string(lp_sol[abs.(lp_sol) .> EPSILON]))
        println("LP Rays are "*string(lp_rays))
        println("====================")
    end

    # [CLEAN] Might actually not need this
    dim = length(lp_sol)
    if length(lp_rays) != dim
        return SCIP.SCIP_DIDNOTFIND 
    end
    
    # STEP 2: Decide Splitting Variable
    vars = get_lp_variables(sepa.scipd)
    split_index = get_first_fractional_index(vars) 
    println(lp_sol[split_index])

    if split_index == -1 
        @warn "No Splitting Variable"
        return SCIP.SCIP_DIDNOTFIND
    end

    intersection_points, parallel_ray = computeIntersectionPoints(lp_sol, split_index, lp_rays)
    
    if DEBUG_PRINT_INTERSECTION_POINTS
        println("====================")
        println("Intersection Points")
        println(intersection_points)
        println("Parallel Ray")
        println(parallel_ray)
        println("====================")
    end

    # STEP 3: Construct and Solve Seperating LP
    lpi = constructSeperatingLP(sepa, dim, lp_sol, intersection_points, parallel_ray) 

    # Solve The LP
    if SCIP.SCIPlpiHasPrimalSolve() == 0
        @warn "LP Solver Does Not Support Primal Solve"
        return SCIP.SCIP_DIDNOTRUN
    end

    SCIP.@SCIP_CALL SCIP.SCIPlpiSolvePrimal(lpi[])
    
    if SCIP.SCIPlpiIsPrimalInfeasible(lpi[]) == 1
        @warn "LP Solver Failed TO Solve Cut Generation LP: LP Infeasible"
        return SCIP.SCIP_DIDNOTFIND
    end
    if SCIP.SCIPlpiIsPrimalUnbounded(lpi[]) == 1
        @warn "LP Solver Failed TO Solve Cut Generation LP: LP Unbounded"
        return SCIP.SCIP_DIDNOTFIND
    end
    if SCIP.SCIPlpiIsOptimal(lpi[]) == 0
        @warn "LP Solver Failed TO Solve Cut Generation LP: LP Not optimzal" 
        return SCIP.SCIP_DIDNOTFIND
    end
    
    obj = Ref{SCIP.SCIP_Real}(0.0)
    reference_sol = zeros(SCIP.SCIP_Real, 2*dim)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(lpi[], obj, pointer(reference_sol),C_NULL, C_NULL, C_NULL)
    reference_sol = reference_sol[1:dim]
    
    #print(sol)
    #println(length(sol))
    b = dot(reference_sol,lp_sol) + 1
    nnonz = 0
    for s in reference_sol
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
    for (index,s) in enumerate(reference_sol)
        if abs(s) >= EPSILON
            cols__[cnt] = cols[index]
            vals__[cnt] = s
            cnt += 1
        end
    end
    
    #println("Largest Coefficient"* string(maximum(abs.(sol))))
    if maximum(abs.(reference_sol)) >= 11000
        println("Coeffients Explotion")
        return SCIP.SCIP_DIDNOTRUN
    end
    SCIP.@SCIP_CALL SCIP.SCIPcreateRowSepa(sepa.scipd, row, sepa.scipd.sepas[sepa], "",nnonz,pointer(cols__),pointer(vals__),b,SCIP.SCIPinfinity(sepa.scipd), true,false, false) 
    #=
    if !isempty(sepa.debug_sol_path)
        println("Debug Solution Available")
        reference_sol = Ref{Ptr{SCIP.SCIP_Sol}}()
        partial = Ref{SCIP.SCIP_Bool}(0)
        error = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPcreateSol(sepa.scipd,reference_sol, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPreadSolFile(sepa.scipd, sepa.debug_sol_path, reference_sol[],SCIP.SCIP_Bool(false),partial, error)
        solution_feasibility = SCIP.SCIPgetRowSolFeasibility(sepa.scipd, row[],reference_sol[])
        println("Solution Feasibility: "*string(solution_feasibility))
        
        if solution_feasibility < -1 
            # In the case of infeasibility
            cols = SCIP.LibSCIP.SCIPgetLPCols(sepa.scipd)
            cols = unsafe_wrap(Vector{Ptr{SCIP.LibSCIP.SCIP_COL}},cols,dim)
            reference_lp = zeros(dim)
            for i = 1:dim
                var = SCIP.LibSCIP.SCIPcolGetVar(cols[i])
                sol = SCIP.SCIPgetSolVal(sepa.scipd, reference_sol[],var)
                name = unsafe_string(SCIP.SCIPvarGetName(var))
                if sol >= EPSILON
                    #println(name* " " * string(sol))
                end
                reference_lp[i] = sol
            end
            #println(reference_lp)
            check_in_corner_polyhedron(lp_sol, lp_rays,reference_lp)
            verify_reference_solution_in_disjunction(sepa,reference_sol,split_index)
        end
    end
    =#
    #SCIP.@SCIP_CALL SCIP.SCIPprintRow(sepa.scipd, row[], C_NULL)
    
    SCIP.@SCIP_CALL SCIP.SCIPaddRow(sepa.scipd,row[],true, infeasible)  
    
    if !isempty(sepa.debug_sol_path)
        println("Debug Solution Available")
        reference_sol = Ref{Ptr{SCIP.SCIP_Sol}}()
        partial = Ref{SCIP.SCIP_Bool}(0)
        error = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPcreateSol(sepa.scipd,reference_sol, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPreadSolFile(sepa.scipd, sepa.debug_sol_path, reference_sol[],SCIP.SCIP_Bool(false),partial, error)
        SCIP.@SCIP_CALL SCIP.SCIPtrySol(sepa.scipd, reference_sol[],1,1,1,1,1,partial)   
        println(partial[] != 0)     
    end
    
    return SCIP.SCIP_SEPARATED
end
