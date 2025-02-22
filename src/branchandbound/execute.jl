using SCIP
using DataStructures
using Dates

function execute_branchandbound(branchandbound::BranchAndBound; log_path = nothing)::Bool
    logger = isnothing(log_path) ? NullLogger() : setup_file_logger(log_path)
    with_logger(logger) do
        _execute_branchandbound(branchandbound)
    end
end

function _execute_branchandbound(branchandbound::BranchAndBound)::Bool
    start_time = time()
    scip = get_scip(branchandbound)

    # Step 1: Initialization
    # We do everything in probing mode
    SCIP.SCIPstartProbing(scip)
    @debug "Starting Branch and Bound"

    # Create root node and put it in node_list
    root = Node(nothing, nothing, 0)
    node_queue_push!(branchandbound, root)

    # Setup loop
    current_node = nothing
    iteration_count = 0

    # Main Branch and Bound Loop
    while !node_queue_empty(branchandbound)
        if branchandbound._time_limit < time() - start_time
            throw(TimeLimitExceededBranchAndBound())
        end
        # Increment iteration count
        iteration_count += 1

        # Step 2: Select next node
        next_node = node_queue_pop!(branchandbound)
        @debug "Node #$(iteration_count): $(next_node) Depth: $(get_depth(next_node))"
        current_node = switch_node!(scip, current_node, next_node)
        prunable = propagate!(scip)
        if prunable
            deactivate!(current_node)
            continue
        end

        # Step 3: Fathom By Bounding
        # Compute dual bound
        lp_feasible = solve_lp_relaxation(scip)
        if !lp_feasible
            @debug "LP is infeasible pruning"
            deactivate!(current_node)
            continue
        end

        # If LP Objective is greater than lower bound then prune
        lp_objective = SCIP.SCIPgetLPObjval(scip)
        @debug "LP Objective: $(lp_objective)"
        if is_GE(lp_objective, get_primal_bound(branchandbound))
            @debug "Node can be pruned by bounding"
            deactivate!(current_node)
            continue
        end

        # Step 4: Fathom By Optimality
        if is_sol_integral(scip)
            @debug "Solution is integral"
            sol = collect_solution(branchandbound)
            deactivate!(current_node)
            # Integral Nodes are leaf
            push_leaf!(branchandbound, current_node)

            # If new solution is better than lower bound then update
            if is_LT(lp_objective, get_primal_bound(branchandbound))
                @debug "New Best Solution Found"
                @debug "Objective: $(lp_objective)"
                set_best_solution(branchandbound, sol)
                set_primal_bound(branchandbound, lp_objective)
            end
            continue
        end

        # Check leaf count limit before Branching (since branching at leaf count)
        max_leaves = get_max_leaves(branchandbound)
        integral_leaves = get_leaf_count(branchandbound)
        node_leaves = get_node_queue_length(branchandbound)
        total_leaves = integral_leaves + node_leaves

        @debug "Total Leaves So Far: $(total_leaves)"
        if max_leaves > 0 && total_leaves + 2 > max_leaves
            # End search
            push_leaf!(branchandbound, current_node)
            break
        end

        # Step 5: Branch
        var = get_branching_variable(branchandbound, scip)
        lpsol = SCIP.SCIPvarGetLPSol(var)
        name = unsafe_string(SCIP.SCIPvarGetName(var))
        # Warning! Some Node Queuing Procedure may alter the SCIP Internal State (e.g. if you use BestFirst) so store LP Solution before hand   
        @debug "Branching on $(name) with LP Solution: $(lpsol)"
        branch(branchandbound, var, UP, ceil(lpsol), current_node)
        @debug "Created Node $(name)>=$(ceil(lpsol))"
        branch(branchandbound, var, DOWN, floor(lpsol), current_node)
        @debug "Created Node $(name)<=$(floor(lpsol))"
    end

    # Collect Leaves
    while !node_queue_empty(branchandbound)
        node = node_queue_pop!(branchandbound)
        push_leaf!(branchandbound, node)
    end

    SCIP.SCIPendProbing(scip)

    # Step 6: Terminate
    # If x is void then ILP is infeasible
    if isnothing(get_best_solution(branchandbound))
        return false
    end

    return true
end

# We use notation from Achterberg Constraint Integer Programming p.46 Algorithm 3.1 
# Modify SCIP to be at the new node. Return new_node 
function switch_node!(
    scip::SCIP.SCIPData, old_node::Union{Nothing,Node}, new_node::Node
)::Node
    # If new node is root then just activate it
    if isroot(new_node)
        activate!(new_node)
        return new_node
    end
    @debug "Switching from $(old_node) to $(new_node)"

    # Step 1: Find common ancestor
    q_hat = new_node
    while !isactive(q_hat)
        q_hat = get_parent(q_hat)
    end
    @debug "Common Ancestor: $(q_hat)"

    # Step 2 undo all changes from q_hat to old_node
    # mark all nodes from old_node to q_hat as inactive
    q = old_node
    while q != q_hat
        @debug "Undoing Action: $(get_action(q))"
        deactivate!(q)
        q = get_parent(q)
    end
    SCIP.SCIPbacktrackProbing(scip, get_depth(q_hat))

    # Step 3: Apply changes from q_hat to new_node
    # Since order needs to be reversed we first get all the action
    actions = []
    q = new_node
    while q != q_hat
        activate!(q)
        push!(actions, get_action(q))
        q = get_parent(q)
    end
    reverse!(actions)

    # Now Apply the actions
    for action in actions
        @debug "Applying Action: $(action)"
        SCIP.SCIPnewProbingNode(scip)
        do_action(scip, action)
    end

    return new_node
end

# Branching Codes
function branch(
    branchandbound::BranchAndBound, var::Ptr{SCIP.SCIP_VAR}, direction::Direction,
    bound::SCIP.SCIP_Real, current_node::Node
)
    action = Action(var, direction, bound)
    new_node = Node(current_node, action, get_depth(current_node) + 1)
    node_queue_push!(branchandbound, new_node)
end