using SCIP
using DataStructures
using AbstractTrees
using D3Trees

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
    var_name = unsafe_string(SCIP.SCIPvarGetName(var))
    if get_direction(action) == DOWN
        print(io, "$(var_name) <= $(get_bound(action))")
    else
        print(io, "$(var_name) >= $(get_bound(action))")
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
    _left::Union{Node,Nothing}
    _right::Union{Node,Nothing}
    """
    Dual bound
    """
    _dual_bound::SCIP.SCIP_Real
end

function Node(action::Union{Action,Nothing}, parent::Union{Node,Nothing}, depth::Int)::Node
    return Node(false, action, parent, depth, nothing, nothing, -Inf)
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

function isroot(node::Node)
    return isnothing(node._action)
end

function get_action(node::Node)
    return node._action
end

function get_parent(node::Node)
    return node._parent
end

function get_depth(node::Node)
    return node._depth
end

function get_left(node::Node)::Union{Node,Nothing}
    return node._left
end

function set_left(node::Node, left::Union{Node,Nothing})
    node._left = left
end

function get_right(node::Node)::Union{Node,Nothing}
    return node._right
end

function set_right(node::Node, right::Union{Node,Nothing})
    node._right = right
end

function get_dual_bound(node::Node)
    return node._dual_bound
end

function set_dual_bound(node::Node, bound::SCIP.SCIP_Real)
    node._dual_bound = bound
end

function Base.show(io::IO, node::Node)
    if isroot(node)
        print(io, "Root")
    else
        print(io, get_action(node))
    end
end

function AbstractTrees.children(node::Node)
    children = []
    if !isnothing(get_left(node))
        push!(children, get_left(node))
    end
    if !isnothing(get_right(node))
        push!(children, get_right(node))
    end
    return children
end

function print_tree(node::Node)
    tree = D3Tree(node)
    display(tree)
end

function prune!(node::Node)
    deactivate!(node)
    if !isroot(node)
        parent = get_parent(node)
        if get_left(parent) == node
            set_left(parent, nothing)
        else
            set_right(parent, nothing)
        end
    end
end