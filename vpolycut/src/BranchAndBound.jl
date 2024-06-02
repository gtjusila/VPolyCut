import SCIP

include("constants.jl")
include("CornerPolyhedron.jl")
include("utils.jl")

function get_point_ray_collection(scip::SCIP.SCIPData, points, ray,fixed; root = false)
    
    if root == false
        #STEP 0: This is a child node, we first need to solve the LP relaxation

        #STEP 0.1: Construct LP Relaxation
        cutoff = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPflushLP(scip)
        SCIP.@SCIP_CALL SCIP.SCIPconstructLP(scip,cutoff)
        if cutoff[] != 0
            @warn "Cutoff during LP Construction"
            return
        end

        #STEP 0.2: Solve LP Relaxation and Cutoff if possible
        error = Ref{SCIP.SCIP_Bool}(C_NULL) 
        cutoff = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPsolveProbingLP(scip, -1, error, cutoff)
        if error[] != 0
            if DEBUG_BRANCH_AND_BOUND
                println("Error During LP Solve")
            end
            return
        end
        if cutoff[] != 0
            if DEBUG_BRANCH_AND_BOUND
                println("====================")
                println("LP Infeasible")
                println("====================")
            end
            return
        end

        # [CLEAN] For Debugging Only
        # SCIP.SCIPwriteLP(scip,joinpath(WRITE_PATH, name*".lp"))

    end

    #Step 1: Select the first somewhat fractional variable
    # [CLEAN] k added temporarily for debugging
    var,sol,k = get_best_branching_candidate(scip; fixed = fixed) 

    #Step 2: Branch
    depth = SCIP.SCIPgetProbingDepth(scip)

    #Step 2.1: Case LEAF
    if depth >= DEPTH_LIMIT || isnothing(var)
        po, ra = get_corner_polyhedron(scip)
        if DEBUG_BRANCH_AND_BOUND
            println("====================")
            println("Leaf Reached. Collecting point and rays")
            println("Point is "*string(po))
            println("Rays are "*string(ra))
            println("====================")
            print_lp_information(scip)
        end
        push!(points,po)
        append!(ray, ra)

        if length(ra) != length(po)
            @warn "Too few ray"
            print(get_dense_tableau_rows(scip))
            SCIP.@SCIP_CALL SCIP.SCIPwriteLP(scip, joinpath(WRITE_PATH,"badray"*string(sepa.called)*".lp"))
        end

        return 
    end

    # [CLEAN] Only For Debugging NOT;A LEAF
    solution = get_lp_solution_vector(scip)
    if DEBUG_BRANCH_AND_BOUND
        println("====================")
        println("Branching")
        println("Current LP Solution: " * string(solution))
        println("Index of branching variable: "*string(k))
        println("Node Depth: "*string(depth)) 
        println("Changing Lower Bound to: "*string(ceil(sol)))
        println("====================")
    end

    #Step 2.2: Go LEFT
    lower = ceil(sol)
    SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(scip)
    SCIP.@SCIP_CALL SCIP.SCIPchgVarLbProbing(scip, var, lower)
    get_point_ray_collection(scip, points, ray, fixed)
    SCIP.@SCIP_CALL SCIP.SCIPbacktrackProbing(scip, depth)
    
    #Step 2.3: GO RIGHT
    
    #[CLEAN] For debugging only
    if DEBUG_BRANCH_AND_BOUND
        println("====================")
        println("Reverted Lower Bound ")
        println("Current LP Solution: " * string(solution))
        println("Index of branching variable: "*string(k))
        println("Node Depth: "*string(depth)) 
        println("Changing Upper Bound to: "*string(floor(sol)))
        println("====================")
    end

    # Actual code starts here
    upper = floor(sol)
    SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(scip)
    SCIP.@SCIP_CALL SCIP.SCIPchgVarUbProbing(scip, var, upper)
    
    # [CLEAN] this is also only temporary for debugging
    if abs(upper) <= 0.5 
        cutoff = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPconstructLP(scip,cutoff)
        SCIP.@SCIP_CALL SCIP.SCIPwriteLP(scip,joinpath(WRITE_PATH, "xlp.lp"))
    end

    # The code continues here
    get_point_ray_collection(scip, points, ray, fixed)
    SCIP.@SCIP_CALL SCIP.SCIPbacktrackProbing(scip, depth)
    
    # Done
    if DEBUG_BRANCH_AND_BOUND
        println("====================")
        println("Node Processing Ended")
        println("Current LP Solution: " * string(solution))
        println("Index of branching variable: "*string(k))
        println("Node Depth: "*string(depth)) 
        println("====================")
    end
end

function get_best_branching_candidate(scip; fixed = [])

    cols = get_lp_columns(scip)
    n = length(cols)
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
    
    #[CLEAN] k added temporarily for debugging
    return var,sol,k
end
