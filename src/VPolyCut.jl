using SCIP
using Random
using LinearAlgebra
using Printf

include("VPolySeparator.jl")
include("debug_utils.jl")
include("CornerPolyhedron.jl")

DEPTH_LIMIT = 2
EPSILON = 1e-6

DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON = true
DEBUG_BRANCH_AND_BOUND = true

WRITE_PATH = joinpath(pwd(),"temp")

function printCurrentLPSolution(sepa::VPolySeparator)
    cols_ = sepa.lp_cols
    n_ = sepa.col_num
    current_solution = zeros(n_)
    for i=1:n_
        var = SCIP.LibSCIP.SCIPcolGetVar(cols_[i])
        current_solution[i] = SCIP.LibSCIP.SCIPvarGetLPSol(var)
    end
    return current_solution
end

function printCut(coef,b)
    str = ""
    for (index,i) in enumerate(coef)
        str = str*" $i"*"x$index"
    end
    str = str*" ≥ $b"
    return str 
end

function printLPtableau(sepa::VPolySeparator)
    n_ =  sepa.col_num
    m_ =  sepa.row_num
    ray_collection_ = Vector{Vector{Int64}}(undef,0) 
    
    # Determine Basic Indices
    basis_indices_ = zeros(Cint,m_)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(sepa.scipd, pointer(basis_indices_))
    println(basis_indices_)
    #Collect Tableau Rows
    tableau_ = Dict()
    for (k, index) in enumerate(basis_indices_)
        if index >= 0
            tableau_[index + 1] = zeros(m_ + n_)
            SCIP.@SCIP_CALL SCIP.LibSCIP.SCIPgetLPBInvRow(sepa.scipd, k-1, pointer(tableau_[index + 1],n_+1), C_NULL, C_NULL)
            SCIP.@SCIP_CALL SCIP.LibSCIP.SCIPgetLPBInvARow(sepa.scipd, k-1, pointer(tableau_[index + 1],n_+1), pointer(tableau_[index + 1]),C_NULL, C_NULL)
            println(tableau_[index + 1])
        end
    end
end

function get_point_ray_collection(sepa::VPolySeparator, points, ray,fixed)
    # If SCIP is not in Probing Mode then start probing model
    if SCIP.SCIPinProbing(sepa.scipd) == 0
        if DEBUG_BRANCH_AND_BOUND
            println("====================")
            println("Starting Probing Mode")
            println("====================")
        end
        SCIP.@SCIP_CALL SCIP.SCIPstartProbing(sepa.scipd)
    else
        #This is a child node and we can collect stuff here
        infeasible = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPconstructLP(sepa.scipd,infeasible)
        if infeasible[] != 0
            @warn "Cutoff during LP Construction"
            return
        end 
        error = Ref{SCIP.SCIP_Bool}(C_NULL) 
        infeasible = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPsolveProbingLP(sepa.scipd, -1,error, infeasible)
        if error[] != 0
            if DEBUG_BRANCH_AND_BOUND
                println("Error During LP Solve")
            end
            return
        end
        if infeasible[] != 0
            if DEBUG_BRANCH_AND_BOUND
                println("====================")
                println("LP Infeasible")
                println("====================")
            end
            return
        end
        name = randstring(3)
        SCIP.SCIPwriteLP(sepa.scipd,joinpath(WRITE_PATH, name*".lp"))
        println("LP written as "*name)
    end
    #Simple Step: Select the first somewhat fractional variable
    cols = sepa.lp_cols
    n = sepa.col_num
    var = nothing
    sol = nothing 
    k = -1
    for i=1:n 
        # Loop through each variable
        var_ = SCIP.SCIPcolGetVar(cols[i])
        sol_ = SCIP.SCIPvarGetLPSol(var_)
        if SCIP.SCIPvarIsIntegral(var_)==1 && (sol_ - floor(sol_) > EPSILON) && !(var_ in fixed)
            #Only Consider The Split if var is integral and the current solution is non integral
            k = i 
            push!(fixed,var_)
            var = var_
            sol = sol_
            break
        end
    end
    depth = SCIP.SCIPgetProbingDepth(sepa.scipd)
    if depth >= DEPTH_LIMIT || isnothing(var)
        po, ra = get_corner_polyhedron(sepa.scipd)
        if DEBUG_BRANCH_AND_BOUND
            println("====================")
            println("Leaf Reached. Collecting point and rays")
            println("Point is "*string(po))
            println("Rays are "*string(ra))
            println("====================")
        end
        push!(points,po)
        for r in ra
            push!(ray,r)
        end
        if length(ra) != length(po)
            @warn "Too few ray"
            printLPtableau(sepa)
            SCIP.@SCIP_CALL SCIP.SCIPwriteLP(sepa.scipd, joinpath(WRITE_PATH,"badray"*string(sepa.called)*".lp"))
        end
        return 
    end
    #Now Upper bound
    lower = ceil(sol)
    upper = floor(sol)
    solution = printCurrentLPSolution(sepa) #[CLEAN] only for debugging
    if DEBUG_BRANCH_AND_BOUND
        println("====================")
        println("Branching")
        println("Current LP Solution: " * string(solution))
        println("Index of branching variable: "*string(k))
        println("Node Depth: "*string(depth)) 
        println("Changing Lower Bound to: "*string(lower))
        println("====================")
    end
    SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(sepa.scipd)
    SCIP.@SCIP_CALL SCIP.SCIPchgVarLbProbing(sepa.scipd, var, lower)
    get_point_ray_collection(sepa, points, ray, fixed)
    SCIP.@SCIP_CALL SCIP.SCIPbacktrackProbing(sepa.scipd, depth)
    if DEBUG_BRANCH_AND_BOUND
        println("====================")
        println("Reverted Lower Bound ")
        println("Current LP Solution: " * string(solution))
        println("Index of branching variable: "*string(k))
        println("Node Depth: "*string(depth)) 
        println("Changing Upper Bound to: "*string(upper))
        println("====================")
    end
    SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(sepa.scipd)
    SCIP.@SCIP_CALL SCIP.SCIPchgVarUbProbing(sepa.scipd, var, upper)
    
    if abs(upper) <= EPSILON
        add_probing_bound_via_row(sepa, var, upper = upper)
    else
        SCIP.@SCIP_CALL SCIP.SCIPchgVarUbProbing(sepa.scipd, var, upper)
    end
    
    get_point_ray_collection(sepa, points, ray, fixed)
    if abs(upper) <= 0.5 
        infeasible = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPconstructLP(sepa.scipd,infeasible)
        SCIP.@SCIP_CALL SCIP.SCIPwriteLP(sepa.scipd,joinpath(WRITE_PATH, "xlp.lp"))
    end
    SCIP.@SCIP_CALL SCIP.SCIPbacktrackProbing(sepa.scipd, depth)
    if DEBUG_BRANCH_AND_BOUND
        println("====================")
        println("Node Processing Ended")
        println("Current LP Solution: " * string(solution))
        println("Index of branching variable: "*string(k))
        println("Node Depth: "*string(depth)) 
        println("====================")
    end
end

function SCIP.exec_lp(sepa::VPolySeparator)
    #Only limit 1 call 
    if sepa.called >=1
        return SCIP.SCIP_DIDNOTRUN
    end
    sepa.called +=1
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
    SCIP.@SCIP_CALL SCIP.SCIPwriteMIP(sepa.scipd, joinpath(WRITE_PATH,name*".mip"),true, false, true)
    # Get Current LP solution 
    lp_sol, lp_rays = get_corner_polyhedron(sepa.scipd)
    if DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON
        println("====================")
        println("Seperator Called")
        println("LP Solution is "*string(lp_sol))
        println("LP Rays are "*string(lp_rays))
        println("====================")
    end
    dim = length(lp_sol)
    if length(lp_rays) != dim
        @warn "Not enough ray"
        return SCIP.SCIP_DIDNOTRUN 
    end
    ray_matrix = stack(lp_rays)
    #Okay So now we get point ray collection
    points_collection = []
    rays_collection = []
    fixed = Set()
    get_point_ray_collection(sepa,points_collection, rays_collection,fixed)
    @info "FInish collecting" 
    @info points_collection 
    @info rays_collection

    #Start Constructing PRLP
    lpi = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(lpi, C_NULL, "Seperating LP", SCIP.SCIP_OBJSEN_MINIMIZE)

    # Add Variables
    SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(lpi[],dim,ones(dim),-SCIP.SCIPinfinity(sepa.scipd)*ones(dim),SCIP.SCIPinfinity(sepa.scipd)*ones(dim),C_NULL, 0,C_NULL, C_NULL, C_NULL )
    for point in points_collection
        lhs = [1.0]
        rhs = [SCIP.SCIPinfinity(sepa.scipd)]
        point = point - lp_sol
        point = ray_matrix \ point
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
    for ray in rays_collection
        lhs = [0.0]
        rhs = [SCIP.SCIPinfinity(sepa.scipd)]
        ray = ray_matrix \ ray
        @info "Ray" ray
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
    println("CGLP Solution is "*string(sol))
    ray_m = ray_matrix'
    sol = ray_m \ sol
    b = dot(sol,lp_sol) + 1
    @info printCut(sol, b)
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
    infeasible = Ref{SCIP.SCIP_Bool}(0)
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

    @warn "Row ADDED"

    return SCIP.SCIP_SEPARATED
end
