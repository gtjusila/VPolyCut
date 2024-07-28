using SCIP
using DataStructures

function execute_branchandbound(branchandbound::BranchAndBound)::Bool
    scip = get_scip(branchandbound)

    # Step 1: Initialization
    # We do everything in probing mode
    SCIP.SCIPstartProbing(scip)

    # Create root node and put it in node_list
    root = Node(nothing, nothing, 0)
    node_queue_enqueue!(branchandbound, root, 0.0)
    # Initial dual bound at root is the solution of LP
    set_dual_bound(root, SCIP.SCIPgetLPObjval(scip))

    # Setup loop
    current_node = nothing
    iteration_count = 1

    # Main Branch and Bound Loop
    while !node_queue_empty(branchandbound)
        if iteration_count > 10000
            return false
        end
        if iteration_count % 100 == 0
            println(
                "Iteration $(iteration_count) with best solution $(get_primal_bound(branchandbound))"
            )
        end

        # Increment iteration count
        iteration_count += 1

        # Step 2: Select next node
        next_node = node_queue_dequeue!(branchandbound)
        current_node = switch_node!(scip, current_node, next_node)
        prunable = propagate!(scip)
        if prunable
            @info "Pruning node line 35"
            prune!(current_node)
            continue
        end

        # Step 3: Fathom By Bounding
        # Compute dual bound
        lp_feasible = solve_lp_relaxation(scip)
        if !lp_feasible
            @info "pruning node line 44"
            prune!(current_node)
            continue
        end

        # If LP Objective is greater than lower bound then prune
        lp_objective = SCIP.SCIPgetLPObjval(scip)
        if SCIP.SCIPisGE(scip, lp_objective, get_primal_bound(branchandbound)) == 1
            #@info "Phatoming node line 52"
            prune!(current_node)
            continue
        end

        # Step 4: Fathom By Optimality
        if is_sol_integral(scip)
            sol = collect_solution(branchandbound)
            deactivate!(current_node)

            # If new solution is better than lower bound then update
            if SCIP.SCIPisLT(scip, lp_objective, get_primal_bound(branchandbound)) == 1
                set_best_solution(branchandbound, sol)
                set_primal_bound(branchandbound, lp_objective)
            end
            continue
        end

        # Step 5: Branch
        # First Fractional branching
        var = get_branching_variable(scip)
        value = SCIP.SCIPvarGetLPSol(var)

        left_action = Action(var, DOWN, floor(value))
        left_prio = get_dual_bound(scip, left_action)
        if SCIP.SCIPisInfinity(scip, left_prio) == 0
            left_node = Node(left_action, current_node, get_depth(current_node) + 1)
            node_queue_enqueue!(branchandbound, left_node, left_prio)
            set_left(current_node, left_node)
        end

        #=
        prio = SCIP.SCIPcalcNodeselPriority(
            scip, var, SCIP.SCIP_BRANCHDIR_UPWARDS, ceil(value)
        )
        =#
        right_action = Action(var, UP, ceil(value))
        right_prio = get_dual_bound(scip, right_action)
        if SCIP.SCIPisInfinity(scip, right_prio) == 0
            right_node = Node(right_action, current_node, get_depth(current_node) + 1)
            node_queue_enqueue!(branchandbound, right_node, right_prio)
            set_right(current_node, right_node)
        end
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

    # Step 1: Find common ancestor
    q_hat = new_node
    while !isactive(q_hat)
        q_hat = get_parent(q_hat)
    end

    # Step 2 undo all changes from q_hat to old_node
    # mark all nodes from old_node to q_hat as inactive
    q = old_node
    while q != q_hat
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
        SCIP.SCIPnewProbingNode(scip)
        do_action(scip, action)
    end

    return new_node
end

# Temporarily execute action 
function get_dual_bound(scip::SCIP.SCIPData, action::Action)::SCIP.SCIP_Real
    depth = SCIP.SCIPgetProbingDepth(scip)
    SCIP.SCIPnewProbingNode(scip)
    do_action(scip, action)

    pruneable = propagate!(scip)
    if pruneable
        @info "Pruning action $(action)"
        return SCIP.SCIPinfinity(scip)
    end

    feasible = solve_lp_relaxation(scip)
    if !feasible
        @info "Infeasible action $(action)"
        return SCIP.SCIPinfinity(scip)
    end

    lpobj = SCIP.SCIPgetLPObjval(scip)
    SCIP.SCIPbacktrackProbing(scip, depth)

    return lpobj
end