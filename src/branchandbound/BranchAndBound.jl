using SCIP
using DataStructures

# Modified Branch And Bound that collect leaves
@kwdef mutable struct BranchAndBound
    _scip::SCIP.SCIPData
    _primal_bound::SCIP.SCIP_Real
    _original_cols::Vector{Ptr{SCIP.SCIP_COL}}
    _best_solution::Union{Nothing,Vector{SCIP.SCIP_Real}}
    _node_queue::NodeQueue
    _max_leaves::Int
    _leaf_nodes::Vector{Node} = Node[]
    _branching_rule::BranchingRule = PseudoCostBranching()
    _time_limit::Float64 = typemax(Float64)
end

function BranchAndBound(
    scip::SCIP.SCIPData;
    max_leaves = -1,
    branching_rule::BranchingRule = PseudoCostBranching(),
    time_limit = typemax(Float64)
)::BranchAndBound
    return BranchAndBound(
        scip,
        BestFirstQueue(; scip = scip);
        max_leaves = max_leaves,
        branching_rule = branching_rule,
        time_limit = time_limit
    )
end

function BranchAndBound(
    scip::SCIP.SCIPData,
    node_queue::NodeQueue;
    max_leaves = -1,
    branching_rule::BranchingRule = PseudoCostBranching(),
    time_limit = typemax(Float64)
)::BranchAndBound
    # Initialize
    n = SCIP.SCIPgetNLPCols(scip)
    original_cols = SCIP.SCIPgetLPCols(scip)
    original_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, original_cols, n)
    return BranchAndBound(;
        _scip = scip,
        _primal_bound = SCIP.SCIPinfinity(scip),
        _original_cols = original_cols,
        _best_solution = nothing,
        _node_queue = node_queue,
        _max_leaves = max_leaves,
        _branching_rule = branching_rule,
        _time_limit = time_limit
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

function get_best_solution(
    branchandbound::BranchAndBound
)::Union{Nothing,Vector{SCIP.SCIP_Real}}
    return branchandbound._best_solution
end

function set_best_solution(
    branchandbound::BranchAndBound, sol::Union{Nothing,Vector{SCIP.SCIP_Real}}
)
    branchandbound._best_solution = sol
end

function node_queue_push!(
    branchandbound::BranchAndBound, node::Node
)
    node_queue_push!(branchandbound._node_queue, node)
end

function node_queue_pop!(branchandbound::BranchAndBound)::Node
    return node_queue_pop!(branchandbound._node_queue)
end

function node_queue_empty(branchandbound::BranchAndBound)::Bool
    return node_queue_empty(branchandbound._node_queue)
end

function get_node_queue_length(branchandbound::BranchAndBound)::Int
    return length(branchandbound._node_queue)
end

function get_max_leaves(branchandbound::BranchAndBound)::Int
    return branchandbound._max_leaves
end

function get_leaf_count(branchandbound::BranchAndBound)::Int
    return length(branchandbound._leaf_nodes)
end

function push_leaf!(branchandbound::BranchAndBound, node::Node)
    push!(branchandbound._leaf_nodes, node)
end

function get_leaves(branchandbound::BranchAndBound)::Vector{Node}
    return branchandbound._leaf_nodes
end

function get_branching_variable(
    branchandbound::BranchAndBound, scip::SCIP.SCIPData
)::Ptr{SCIP.SCIP_VAR}
    return get_branching_variable(branchandbound._branching_rule, scip)
end