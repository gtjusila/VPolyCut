using JSON
function get_point_ray_collection(
    sepa::VPCSeparator
)
    scip = sepa.scipd
    disjunction = sepa.disjunction
    projection = sepa.projection
    @debug "Collecting Points and Rays from disjunction"

    # We only need to store projected points and ray
    sepa.point_ray_collection = PointRayCollection(scip; projection=projection)

    # A bin for storing disjunction information, will be pushed to the disjunction information function
    disjunctive_terms = Vector{Dict{String,Tuple{SCIP.SCIP_Real,SCIP.SCIP_Real}}}()

    # Iterate through every disjuncive terms
    for (i, term) in enumerate(disjunction)
        @debug "Collecting Point and Ray from term $i"
        path = get_path(term) # Get all actions from root to leaf
        get_disjunctive_term_information(
            sepa, path, disjunctive_terms
        )
    end

    # Write the disjunction that will be used for cut generation 
    disjunctive_terms_path = joinpath(sepa.parameters.log_directory, "disjunction.json")
    # Write the JSON data to the file
    open(disjunctive_terms_path, "w") do io
        JSON.print(io, disjunctive_terms, 4)
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
    path::Vector{Node},
    disjunctive_terms::Vector{Dict{String,Tuple{SCIP.SCIP_Real,SCIP.SCIP_Real}}}
)
    # Some alias 
    root_tableau = sepa.complemented_tableau
    scip = sepa.scipd

    # Enter Probing mode
    SCIP.SCIPstartProbing(scip)
    # We need to store the affected variable for storing disjunction information later
    affected_vars = Vector{Ptr{SCIP.SCIP_Var}}()

    # Get To Node while keeping track of affected variables
    for node in path
        if isroot(node)
            continue
        end
        action = get_action(node)
        var = get_var(action)
        push!(affected_vars, var)
        do_action(scip, action)
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

    # Record disjunction Information
    # We store disjunction as a dictionary where the key is variable name and the value is a tuple of lower bound and upper bound
    disjunction_information = Dict{String,Tuple{SCIP.SCIP_Real,SCIP.SCIP_Real}}()
    # Go trough Affected variable
    for var in affected_vars
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)
        disjunction_information[unsafe_string(SCIP.SCIPvarGetName(var))] = (lb, ub)
    end
    # Push to the bin
    push!(disjunctive_terms, disjunction_information)

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
    orig_objective_value = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    basic_solution = get_lp_sol(corner_polyhedron)

    # Leave Probing mode
    SCIP.SCIPendProbing(scip)

    # Add rays to point ray collection
    add_point(
        sepa.point_ray_collection, basic_solution, objective_value, orig_objective_value
    )
    for ray in get_lp_rays(corner_polyhedron)
        add_ray(sepa.point_ray_collection, ray)
    end

    @debug "Disjunctive term is feasible. Points and rays added. Objective value is $objective_value"
end