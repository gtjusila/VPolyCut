using DataStructures
using Random
abstract type NodeQueue end

function node_queue_push!(node_queue::T, node::Node) where {T<:NodeQueue}
    error("Method node_queue_push! not implemented for type $(T)")
end

function node_queue_pop!(node_queue::T)::Node where {T<:NodeQueue}
    error("Method node_queue_pop! not implemented for type $(T)")
end

function node_queue_empty(node_queue::T)::Bool where {T<:NodeQueue}
    error("Method node_queue_empty not implemented for type $(T)")
end

function Base.length(node_queue::T)::Int where {T<:NodeQueue}
    error("Method length not implemented for type $(T)")
end

@kwdef mutable struct RandomNodeQueue <: NodeQueue
    _queue::PriorityQueue{Node,Float64} = PriorityQueue()
end

function node_queue_push!(node_queue::RandomNodeQueue, node::Node)
    enqueue!(node_queue._queue, node, rand())
end

function node_queue_pop!(node_queue::RandomNodeQueue)::Node
    return dequeue!(node_queue._queue)
end

function node_queue_empty(node_queue::RandomNodeQueue)::Bool
    return isempty(node_queue._queue)
end

function Base.length(node_queue::RandomNodeQueue)::Int
    return length(node_queue._queue)
end

@kwdef mutable struct BFSQueue <: NodeQueue
    _queue::Queue{Node} = Queue{Node}()
end

function node_queue_push!(node_queue::BFSQueue, node::Node)
    enqueue!(node_queue._queue, node)
end

function node_queue_pop!(node_queue::BFSQueue)::Node
    return dequeue!(node_queue._queue)
end

function node_queue_empty(node_queue::BFSQueue)::Bool
    return isempty(node_queue._queue)
end

function Base.length(node_queue::BFSQueue)::Int
    return length(node_queue._queue)
end

@kwdef mutable struct DFSQueue <: NodeQueue
    _queue::Stack{Node} = Stack{Node}()
end

function node_queue_push!(node_queue::DFSQueue, node::Node)
    push!(node_queue._queue, node)
end

function node_queue_pop!(node_queue::DFSQueue)::Node
    return pop!(node_queue._queue)
end

function node_queue_empty(node_queue::DFSQueue)::Bool
    return isempty(node_queue._queue)
end

function Base.length(node_queue::DFSQueue)::Int
    return length(node_queue._queue)
end

mutable struct BestFirstQueue <: NodeQueue
    _queue::PriorityQueue{Node,Float64}
    _scip::SCIP.SCIPData
end

function BestFirstQueue(; scip::SCIP.SCIPData)
    return BestFirstQueue(PriorityQueue(), scip)
end

function get_scip(node_queue::BestFirstQueue)
    return node_queue._scip
end

function node_queue_push!(node_queue::BestFirstQueue, node::Node)
    # Exception for Root Node 
    if isroot(node)
        enqueue!(node_queue._queue, node, 0)
        return nothing
    end

    scip = get_scip(node_queue)
    bound = get_dual_bound(scip, get_action(node))
    if !is_infinity(bound)
        enqueue!(node_queue._queue, node, bound)
    end
end

function node_queue_pop!(node_queue::BestFirstQueue)::Node
    return dequeue!(node_queue._queue)
end

function node_queue_empty(node_queue::BestFirstQueue)::Bool
    return isempty(node_queue._queue)
end

function Base.length(node_queue::BestFirstQueue)::Int
    return length(node_queue._queue)
end