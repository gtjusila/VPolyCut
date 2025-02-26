using SCIP
using DataStructures

@enum Direction begin
    DOWN
    UP
end

mutable struct Action
    """ 
    The variable affected 
    """
    _var::Ptr{SCIP.SCIP_VAR}
    """
    Does the node bound the variable up or down
    """
    _direction::Direction
    """
    The value of the bound
    """
    _bound::SCIP.SCIP_Real
end

function get_bound(action::Action)
    return action._bound
end

function get_direction(action::Action)
    return action._direction
end

function get_var(action::Action)
    return action._var
end

function do_action(scip::SCIP.SCIPData, action::Action)
    var = get_var(action)
    direction = get_direction(action)
    bound = get_bound(action)
    if direction == UP
        SCIP.SCIPchgVarLbProbing(scip, var, bound)
    else
        SCIP.SCIPchgVarUbProbing(scip, var, bound)
    end
end

function Base.show(io::IO, action::Action)
    var = get_var(action)
    direction = get_direction(action)
    bound = get_bound(action)
    name = unsafe_string(SCIP.SCIPvarGetName(var))
    if direction == UP
        print(io, "$(name) >= $bound")
    else
        print(io, "$(name) <= $bound")
    end
end

mutable struct Node
    """
    is node active
    """
    _active::Bool
    """
    Node action
    """
    _action::Union{Action,Nothing}
    """
    Parent is either a node or nothing if root
    """
    _parent::Union{Node,Nothing}
    """
    Extra Information To Fasten Backtracking
    """
    _depth::Int
    """
    Dual bound
    """
    _dual_bound::SCIP.SCIP_Real
end

function Node(parent::Union{Node,Nothing}, action::Union{Action,Nothing}, depth::Int)
    return Node(false, action, parent, depth, typemax(SCIP.SCIP_Real))
end

function isactive(node::Node)
    return node._active
end

function activate!(node::Node)
    node._active = true
end

function deactivate!(node::Node)
    node._active = false
end

function get_action(node::Node)
    return node._action
end

function get_parent(node::Node)
    return node._parent
end

function isroot(node::Node)
    return isnothing(node._action)
end

function get_depth(node::Node)
    return node._depth
end

function get_path(node::Node)::Vector{Node}
    path = []
    current = node
    while !isroot(current)
        push!(path, current)
        current = get_parent(current)
    end
    push!(path, current)
    return collect(reverse(path))
end

function print_path(node::Node)
    path = get_path(node)
    for node in path
        print(node, ", ")
    end
    print("\n")
end

function Base.show(io::IO, node::Node)
    if isroot(node)
        print(io, "Root")
    else
        action = get_action(node)
        print(io, action)
    end
end