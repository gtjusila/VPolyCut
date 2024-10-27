function get_point_ray_collection(
    sepa::VPCSeparator
)
    scip = sepa.scipd
    disjunction = sepa.disjunction
    projection = sepa.projection
    @debug "Collecting Points and Rays from disjunction"

    # We only need to store projected points and ray
    sepa.point_ray_collection = PointRayCollection(scip; projection=projection)
    for (i, term) in enumerate(disjunction)
        @debug "Collecting Point and Ray from term $i"
        path = get_path(term) # Get all actions from root to leaf
        get_disjunctive_term_information(
            sepa, path
        )
    end

    # Clean Up Duplicate Ray
    @debug "Removing duplicate rays"
    remove_duplicate_rays(sepa.point_ray_collection)
    @debug "Points and Rays Collected"
end

"""
The the point ray collection and solution value are collected from the disjunction 
The point ray collection will be given in the complemented space therefore the 
original root_tableau is required. Point and rays found are added to the sepa point ray collection directly
"""
function get_disjunctive_term_information(
    sepa::VPCSeparator,
    path::Vector{Node}
)
    # Some alias 
    root_tableau = sepa.complemented_tableau
    scip = sepa.scipd

    # Enter Probing mode
    SCIP.SCIPstartProbing(scip)

    # Get To Node
    for node in path
        if isroot(node)
            continue
        end
        do_action(scip, get_action(node))
    end

    # Propegate
    prunable = propagate!(scip)
    if prunable
        @debug "Disjunctive term is detected infeasible during propagation"
        SCIP.SCIPendProbing(scip)
        return nothing, [], 0
    end

    # Solve LP
    lp_feasible = solve_lp_relaxation(scip)
    if !lp_feasible
        @debug "Disjunctive term is detected infeasible during LP solve"
        SCIP.SCIPendProbing(scip)
        return nothing, [], 0
    end

    # Get Optimal Tableau
    tableau = construct_tableau(scip)

    # Complement the same columns as the root tableau
    complemented_columns = get_complemented_columns(root_tableau)
    for i in complemented_columns
        var = get_var_from_column(tableau, i)
        complement_column(tableau, var)
    end

    corner_polyhedron = construct_corner_polyhedron(tableau)
    objective_value = SCIP.SCIPgetLPObjval(scip)
    basic_solution = get_lp_sol(corner_polyhedron)

    # Leave Probing mode
    SCIP.SCIPendProbing(scip)

    # Add rays to point ray collection
    add_point(sepa.point_ray_collection, basic_solution, objective_value)
    for ray in get_lp_rays(corner_polyhedron)
        add_ray(sepa.point_ray_collection, ray)
    end

    @debug "Disjunctive term is feasible. Points and rays added. Objective value is $objective_value"
end