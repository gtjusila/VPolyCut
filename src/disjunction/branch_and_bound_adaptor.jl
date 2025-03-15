using SCIP

"""
create a partial branch and bound tree and return a disjunction where each term is a leaf node of the tree
"""
function get_disjunction_by_branchandbound(
    sepa::VPCSeparator;
)
    scipd = sepa.shared_data.scipd
    n_leaves = sepa.parameters.n_leaves
    log_path = joinpath(sepa.parameters.log_directory, "branch_and_bound.log")
    time_limit = 0.5 * sepa.parameters.time_limit

    # First we generate a branch and bound tree with n_leaves
    branchandbound = BranchAndBound(
        scipd; max_leaves = n_leaves, time_limit = time_limit
    )
    starting_lp_iter_count = SCIP.SCIPgetNLPIterations(scipd)
    execute_branchandbound(branchandbound; log_path = log_path)
    end_lp_iter_count = SCIP.SCIPgetNLPIterations(scipd)

    sepa.statistics.point_ray_collection_lp_iterations =
        end_lp_iter_count - starting_lp_iter_count

    # Collect the leaves
    leaves = get_leaves(branchandbound)

    if isempty(leaves)
        throw(TreeHasNoLeaf())
    end

    # Prepare the disjunction object
    disjunction = Disjunction()

    # For each leaf we create a disjunctive term
    for leaf in leaves
        term = DisjunctiveTerm()
        path = get_path(leaf)

        # We get all the bound changes (action) that is applied from the root to the leaf
        for node in path
            if isroot(node)
                continue
            end
            node_action = get_action(node)
            action_direction = get_direction(node_action)
            action_var = get_var(node_action)
            action_bound = get_bound(node_action)
            if action_direction == UP
                set_lower_bound!(term, action_var, action_bound)
            elseif action_direction == DOWN
                set_upper_bound!(term, action_var, action_bound)
            end
        end

        add_term!(disjunction, term)
    end
    return disjunction
end