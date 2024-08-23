abstract type BranchingRule end

function get_branching_variable(
    branching::T, scip::SCIP.SCIPData
) where {T<:BranchingRule}
    error("BranchingRule not implemented for type $(T)")
end

struct PseudoCostBranching <: BranchingRule end;

function get_branching_variable(
    branching::PseudoCostBranching, scip::SCIP.SCIPData
)::Ptr{SCIP.SCIP_VAR}
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

struct FirstFractionalBranching <: BranchingRule end;

function get_branching_variable(
    branching::FirstFractionalBranching, scip::SCIP.SCIPData
)::Ptr{SCIP.SCIP_VAR}
    # Step 1: Get Branching Candidate
    lpcands = Ref{Ptr{Ptr{SCIP.SCIP_VAR}}}(C_NULL)
    nlpcands = Ref{Cint}(0)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBranchCands(
        scip, lpcands, C_NULL, C_NULL, nlpcands, C_NULL, C_NULL
    )
    lpcands = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, lpcands[], nlpcands[])

    return lpcands[1]
end

struct PriorityBranching <: BranchingRule
    _priority_var::Ptr{SCIP.SCIP_VAR}
end;

function get_priority_var(branching::PriorityBranching)::Ptr{SCIP.SCIP_VAR}
    return branching._priority_var
end

function get_branching_variable(
    branching::PriorityBranching, scip::SCIP.SCIPData
)::Ptr{SCIP.SCIP_VAR}
    prio_var = get_priority_var(branching)
    if SCIP.SCIPvarIsInLP(prio_var) == 1 &&
        SCIP.SCIPvarIsIntegral(prio_var) == 1 &&
        !is_integral(scip, SCIP.SCIPvarGetLPSol(prio_var))
        return prio_var
    end
    return get_branching_variable(PseudoCostBranching(), scip)
end