function solve_separation_subproblems(sepa::VPCSeparator)
    # Setup aliases
    scip = sepa.scipd

    # Get current LP Solution projected onto the non-basic variables.
    # This will be the origin point 
    current_lp_solution = get_solution_vector(sepa.complemented_tableau)
    projected_current_lp_solution = project(sepa.projection, current_lp_solution)
    problem_dimension = length(projected_current_lp_solution)

    # Start building the model
    @debug "Using simplex strategy $(sepa.parameters.lp_solving_method) for HiGHS"
    @debug "Zeroing heuristic is $(sepa.parameters.zeroing_heuristic)"
    separating_lp = Model(HiGHS.Optimizer)
    JuMP.set_optimizer_attribute(separating_lp, "output_flag", false)
    JuMP.set_optimizer_attribute(
        separating_lp, "simplex_strategy", sepa.parameters.lp_solving_method
    )

    # First, get the point ray collection
    # The points should be translated to the non-basic space (they are already projected)
    projected_point_collection = get_points(sepa.point_ray_collection)
    point_collection_in_nonbasic_space = map(projected_point_collection) do corner_point
        return CornerPoint(
            get_point(corner_point) - projected_current_lp_solution,
            get_objective_value(corner_point),
            get_orig_objective_value(corner_point)
        )
    end

    # The rays have already been projected to the non-basic space (rays does not need to be projected)
    ray_collection = get_rays(sepa.point_ray_collection)

    # Construct PRLP0
    # Variables
    @variable(separating_lp, x[1:problem_dimension])

    # Constraints for Points
    for point in point_collection_in_nonbasic_space
        point_coordinates = get_point(point)
        @constraint(
            separating_lp,
            sum(x[i] * point_coordinates[i] for i in 1:problem_dimension) >= 1
        )
    end

    # Constraints for Rays
    for ray in ray_collection
        ray_coefficients = get_coefficients(ray)
        @constraint(
            separating_lp, sum(x[i] * ray_coefficients[i] for i in 1:problem_dimension) >= 0
        )
    end

    # Compute the disjunctive lower bound and the corresponding point p_star
    p_star = argmin(point -> get_objective_value(point), point_collection_in_nonbasic_space)
    sepa.disjunctive_lower_bound = get_orig_objective_value(p_star)
    @debug "Disjunctive Lower Bound is: $(sepa.disjunctive_lower_bound) and current LP Objective is: $(sepa.lp_obj)"
    if !is_GT(scip, sepa.disjunctive_lower_bound, sepa.lp_obj)
        @debug "Disjunctive Lower Bound is not greater than current LP Objective"
        throw(FailedDisjunctiveLowerBoundTest())
    end

    # Execute Cbar Test For Numerical Stability Check
    sepa_tableau_variables = map(1:get_nvars(sepa.complemented_tableau)) do i
        return get_var_from_column(sepa.complemented_tableau, i)
    end
    cbar = map(sepa_tableau_variables) do var
        return get_obj(var)
    end
    projected_cbar = project(sepa.projection, cbar)
    control_solution_coordinates = projected_cbar / dot(projected_cbar, get_point(p_star))
    control_solution_dict = Dict{typeof(x[1]),Float64}(
        x[i] => control_solution_coordinates[i] for i in 1:problem_dimension
    )
    control_solution_feasibility_report = primal_feasibility_report(
        separating_lp, control_solution_dict; atol=1e-7
    )
    if (length(control_solution_feasibility_report) > 0)
        throw(FailedCbarTest())
    end

    # check if the LP in solver feasible by optimizing using all 0s objective
    # Time limit is 300 s
    set_time_limit_sec(separating_lp, 300.0)
    @objective(separating_lp, Min, 0)
    @debug "Starting Check Feasibility"
    optimize!(separating_lp)
    push!(sepa.prlp_solves, get_solve_stat(separating_lp, "feasibility"))
    if primal_status(separating_lp) != MOI.FEASIBLE_POINT
        throw(FailedToProvePRLPFeasibility())
    end
    @debug "Feasibility Check Passed"

    # Now optimize with all 1s objective
    # Time limit is 30 s
    set_time_limit_sec(separating_lp, 30.0)
    @objective(separating_lp, Min, sum(x))
    optimize!(separating_lp)
    push!(sepa.prlp_solves, get_solve_stat(separating_lp, "all_ones"))
    if is_solved_and_feasible(separating_lp)
        cut = get_cut_from_separating_solution(sepa, value.(x))
        push!(sepa.cutpool, cut)
    else
        @warn "All ones objective is unbounded"
    end

    # Now Optimize with p_star as objective
    @debug "Starting P* Optimization"
    @objective(separating_lp, Min, sum(x[i] * p_star[i] for i in 1:problem_dimension))
    optimize!(separating_lp)
    push!(sepa.prlp_solves, get_solve_stat(separating_lp, "p_star"))
    if is_solved_and_feasible(separating_lp)
        cut = get_cut_from_separating_solution(sepa, value.(x))
        push!(sepa.cutpool, cut)
    else
        # Cannot proceed if P* is unbounded
        throw(PStarUnbounded())
    end
    @debug "Finished P* Optimization"

    # Only proceed if the cut is tight at p_star
    if !is_EQ(scip, objective_value(separating_lp), 1.0)
        throw(PStarNotTight())
    end

    # Now transition to PRLP=, add constraint a_bar^Tp_star == 1 
    a_bar = value.(x)
    @constraint(separating_lp, sum(x[i] * p_star[i] for i in 1:problem_dimension) == 1)

    r_bar = filter(ray -> !is_zero(scip, dot(a_bar, ray)), ray_collection)
    sort!(r_bar; by=ray -> abs(get_obj(get_generating_variable(ray))))
    objective_tried = 0

    start_time = time()
    for ray in r_bar
        @objective(
            separating_lp, Min, sum(x[i] * ray[i] for i in 1:problem_dimension)
        )
        optimize!(separating_lp)

        objective_tried += 1
        push!(sepa.prlp_solves, get_solve_stat(separating_lp, "prlp=$(objective_tried)"))
        @info "Objective tried: $(objective_tried). Limit is $(2 * sepa.parameters.cut_limit)"

        if is_solved_and_feasible(separating_lp)
            cut = get_cut_from_separating_solution(sepa, value.(x))
            push!(sepa.cutpool, cut)
            @debug "Found a cut. Cutpool size: $(length(sepa.cutpool)) cuts. Cutlimit: $(sepa.parameters.cut_limit)"
        else
            @debug "Failed optimizing over objective."
        end

        # Termination Condition
        if objective_tried > 2 * sepa.parameters.cut_limit
            break
        end
        if length(sepa.cutpool) >= sepa.parameters.cut_limit
            break
        end
        if time() - start_time > 1800
            break
        end
    end

    if length(sepa.cutpool) > 0
        sepa.separated = true
    end

    add_all_cuts!(sepa.cutpool, sepa)
end

function zeroing_heuristic(
    vector::Vector{SCIP.SCIP_Real}
)::Union{Vector{SCIP.SCIP_Real},Int}
    zeroed = 0
    vector = map(vector) do x
        if is_zero(scip, x)
            if x != 0.0
                zeroed += 1
            end
            return 0.0
        else
            return x
        end
    end
    return vector, zeroed
end

function get_cut_from_separating_solution(
    sepa::VPCSeparator,
    separating_solution::Vector{SCIP.SCIP_Real}
)::Cut
    scip = sepa.scipd
    tableau = sepa.complemented_tableau
    lp_solution = get_solution_vector(tableau)
    lp_solution = get_uncomplemented_vector(lp_solution, tableau)

    separating_solution = undo_projection(sepa.projection, separating_solution)
    separating_solution = get_uncomplemented_vector(separating_solution, tableau)
    b = dot(separating_solution, lp_solution) + 1

    cut_vector, b = convert_standard_inequality_to_general(
        scip, tableau, separating_solution, b
    )

    # We normalize the cut to the form ax <= b
    return Cut(-cut_vector, -b)
end

function get_solve_stat(model::JuMP.AbstractModel, obj_name::String)
    return Dict(
        "solve_time" => solve_time(model),
        "simplex_iterations" => simplex_iterations(model),
        "primal_status" => primal_status(model),
        "dual_status" => dual_status(model),
        "termination_status" => termination_status(model),
        "objective_name" => obj_name
    )
end