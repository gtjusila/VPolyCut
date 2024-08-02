# Miscelaneous function that connect SCIP with Julia
using SCIP
using DataStructures

"""
Propegate bound change within SCIP. Return whether or not bound is pruneable
"""
function propagate!(scip::SCIP.SCIPData)::Bool
    cutoff = Ref{SCIP.SCIP_Bool}(0)
    SCIP.SCIPpropagateProbing(scip, -1, cutoff, C_NULL)
    return Bool(cutoff[])
end

"""
 Solve LP relaxation of current node and return whether or not LP is feasible 
"""
function solve_lp_relaxation(scip::SCIP.SCIPData)::Bool
    if SCIP.SCIPisLPConstructed(scip) == 0
        lp_infeasible = Ref{SCIP.SCIP_Bool}(0)
        SCIP.@SCIP_CALL SCIP.SCIPconstructLP(scip, lp_infeasible)
        lp_infeasible = Bool(lp_infeasible[])
        # If node can be cutoff then prune
        if lp_infeasible
            return false
        end
    end

    # Solve LP
    lp_infeasible = Ref{SCIP.SCIP_Bool}(0)
    lperror = Ref{SCIP.SCIP_Bool}(0)
    SCIP.@SCIP_CALL SCIP.SCIPsolveProbingLP(scip, -1, lperror, lp_infeasible)
    lp_infeasible = Bool(lp_infeasible[])
    lperror = Bool(lperror[])

    if lperror
        error("Error in solving LP. You need to take a look at this")
        exit(1)
    end

    if lp_infeasible
        return false
    end

    return true
end

function is_sol_integral(scip::SCIP.SCIPData)::Bool
    return SCIP.SCIPgetNLPBranchCands(scip) == 0
end

function collect_solution(
    branchandbound::BranchAndBound
)::Vector{Float64}
    return [
        sol_from_column(col) for col in get_original_cols(branchandbound)
    ]
end

function sol_from_column(col::Ptr{SCIP.SCIP_Col})::SCIP.SCIP_Real
    var = SCIP.SCIPcolGetVar(col)
    return SCIP.SCIPvarGetLPSol(var)
end

function get_branching_variable(scip::SCIP.SCIPData)::Ptr{SCIP.SCIP_VAR}
    # Step 1: Get Branching Candidate
    lpcands = Ref{Ptr{Ptr{SCIP.SCIP_VAR}}}(C_NULL)
    nlpcands = Ref{Cint}(0)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBranchCands(
        scip, lpcands, C_NULL, C_NULL, nlpcands, C_NULL, C_NULL
    )
    lpcands = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, lpcands[], nlpcands[])

    # Step 2: Get Pseudocost Score
    priority_queue = PriorityQueue{Ptr{SCIP.SCIP_VAR},Float64}()
    for var in lpcands
        score = SCIP.SCIPgetVarPseudocostScore(scip, var, SCIP.SCIPvarGetLPSol(var))
        enqueue!(priority_queue, var, -score)
    end
    var = dequeue!(priority_queue)

    return var
end

# Temporarily execute action 
function get_dual_bound(scip::SCIP.SCIPData, action::Action)::SCIP.SCIP_Real
    depth = SCIP.SCIPgetProbingDepth(scip)
    SCIP.SCIPnewProbingNode(scip)
    do_action(scip, action)

    pruneable = propagate!(scip)
    if pruneable
        SCIP.SCIPbacktrackProbing(scip, depth)
        return SCIP.SCIPinfinity(scip)
    end

    feasible = solve_lp_relaxation(scip)
    if !feasible
        SCIP.SCIPbacktrackProbing(scip, depth)
        return SCIP.SCIPinfinity(scip)
    end

    lpobj = SCIP.SCIPgetLPObjval(scip)
    SCIP.SCIPbacktrackProbing(scip, depth)
    return lpobj
end