using SCIP
using DataStructures

mutable struct BranchAndBound
    _scip::SCIP.SCIPData
    _primal_bound::SCIP.SCIP_Real
    _original_cols::Vector{Ptr{SCIP.SCIP_COL}}
    _best_solution::Union{Nothing,Point}
    _node_queue::NodeQueue
end

function BranchAndBound(scip::SCIP.SCIPData)
    return BranchAndBound(
        scip,
        BestFirstQueue(; scip=scip)
    )
end
function BranchAndBound(scip::SCIP.SCIPData, node_queue::NodeQueue)::BranchAndBound
    # Initialize
    n = SCIP.SCIPgetNLPCols(scip)
    original_cols = SCIP.SCIPgetLPCols(scip)
    original_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, original_cols, n)

    return BranchAndBound(
        scip,
        SCIP.SCIPinfinity(scip),
        original_cols,
        nothing,
        node_queue
    )
end

function get_scip(branchandbound::BranchAndBound)::SCIP.SCIPData
    return branchandbound._scip
end

function get_primal_bound(branchandbound::BranchAndBound)::SCIP.SCIP_Real
    return branchandbound._primal_bound
end

function set_primal_bound(branchandbound::BranchAndBound, bound::SCIP.SCIP_Real)
    branchandbound._primal_bound = bound
end

function get_original_cols(branchandbound::BranchAndBound)::Vector{Ptr{SCIP.SCIP_COL}}
    return branchandbound._original_cols
end

function get_best_solution(branchandbound::BranchAndBound)::Union{Nothing,Point}
    return branchandbound._best_solution
end

function set_best_solution(branchandbound::BranchAndBound, sol::Union{Nothing,Point})
    branchandbound._best_solution = sol
end

function node_queue_push(
    branchandbound::BranchAndBound, node::Node
)
    node_queue_push!(branchandbound._node_queue, node)
end

function node_queue_pop(branchandbound::BranchAndBound)::Node
    return node_queue_pop!(branchandbound._node_queue)
end

function node_queue_empty(branchandbound::BranchAndBound)::Bool
    return node_queue_empty(branchandbound._node_queue)
end