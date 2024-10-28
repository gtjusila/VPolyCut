function solve_separation_subproblems(sepa::VPCSeparator)
    scip = sepa.scipd

    # Get current LP Solution projected onto the non-basic variables.
    # This will be the origin point 
    lp_solution = get_solution_vector(sepa.complemented_tableau)
    lp_solution = project(sepa.projection, lp_solution)
    dimension = length(lp_solution)

    #Start building the model
    separating_lp = Model(HiGHS.Optimizer)
    JuMP.set_optimizer_attribute(separating_lp, "output_flag", false)
    JuMP.set_optimizer_attribute(separating_lp, "simplex_strategy", 4)
    @debug "Zeroing heuristic is $(sepa.parameters.zeroing_heuristic)"

    # First, translate the points so that the lp_solution is at the origin 
    points = get_points(sepa.point_ray_collection)
    translated_points = map(points) do point
        return CornerPoint(
            get_point(point) - lp_solution,
            get_objective_value(point),
            get_orig_objective_value(point)
        )
    end
    rays = get_rays(sepa.point_ray_collection)
    zeroed = 0

    # Construct PRLP0

    # Variables
    @variable(separating_lp, x[1:dimension])

    # Constraints for Points
    for point in translated_points
        vertex = get_point(point)
        if sepa.parameters.zeroing_heuristic
            vertex, zero_added = zeroing_heuristic(vertex)
            zeroed += zero_added
        end
        @constraint(separating_lp, sum(x[i] * vertex[i] for i in 1:dimension) >= 1)
    end

    # Constraints for Rays
    for ray in rays
        coefficient = get_coefficients(ray)
        if sepa.parameters.zeroing_heuristic
            coefficient, zero_added = zeroing_heuristic(coefficient)
            zeroed += zero_added
        end
        @constraint(separating_lp, sum(x[i] * coefficient[i] for i in 1:dimension) >= 0)
    end

    zeroed = 0
    vars = map(1:get_nvars(sepa.complemented_tableau)) do i
        return get_var_from_column(sepa.complemented_tableau, i)
    end
    c_bar = map(vars) do var
        return get_obj(var)
    end
    lp_sol = map(vars) do var
        return get_sol(var)
    end

    c_bar = project(sepa.projection, c_bar)
    p_star = argmin(point -> get_objective_value(point), translated_points)
    p_star = get_point(p_star)
    frac = dot(c_bar, p_star)
    check = c_bar / frac
    checkvar = Dict{typeof(x[1]),Float64}(x[i] => check[i] for i in 1:dimension)

    # Open the file in write mode
    #=
    open("primal_feasibility_report_filtered.txt", "w") do file
        # Iterate over the dictionary
        p = primal_feasibility_report(separating_lp, checkvar; atol=1e-7)
        if (length(p) > 0)
            write(file, "Complemented\n")
            complemented_columns = get_complemented_columns(sepa.complemented_tableau)
            for i in get_complemented_columns(sepa.complemented_tableau)
                write(file, "$i \n")
            end
            write(file, "Check Variables\n")
            for i in 1:dimension
                # Only write entries where the absolute value of the value is greater than atol
                p = sepa.projection.map_projected_to_original[i]
                complemented = p in complemented_columns
                original = get_var_from_column(sepa.complemented_tableau, p)
                write(
                    file,
                    string(
                        "x[$i] => ",
                        checkvar[x[i]],
                        " is complemented? ",
                        complemented,
                        " original $p",
                        "ub $(get_ub(original)) lb $(get_lb(original)) red $(get_obj(original)) type $(typeof(original)) status $(get_basis_status(original))\n"
                    )
                )
            end
            write(file, "Primal Feasibility Report\n")
            for (key, value) in
                primal_feasibility_report(separating_lp, checkvar; atol=1e-7)
                # Only write entries where the absolute value of the value is greater than atol
                write(file, string(key, " => ", value, "\n"))
            end
        else
            @debug "Check solution is feasible"
        end
    end
    =#
    @debug "Zeroed $(zeroed) entries"

    # Check if disjunctive lower bound is worse than optimal LP objective
    p_star = argmin(point -> get_objective_value(point), translated_points)
    disjunctive_lower_bound = get_objective_value(p_star)
    sepa.disjunctive_lower_bound = get_orig_objective_value(p_star)
    @debug "Disjunctive Lower Bound is: $(disjunctive_lower_bound) and LP Objective is: $(sepa.lp_obj)"
    if !is_GT(scip, disjunctive_lower_bound, sepa.lp_obj)
        throw(FailedDisjunctiveLowerBoundTest())
    end

    # check if the LP is feasible by optimizing using all 0s objective
    @objective(separating_lp, Min, 0)
    write_to_file(separating_lp, "separating_lp.lp")
    @debug "Starting Check Feasibility"
    set_time_limit_sec(separating_lp, 1000.0)
    optimize!(separating_lp)
    if primal_status(separating_lp) != MOI.FEASIBLE_POINT
        throw(FailedToProvePRLPFeasibility())
    end
    @debug "Feasibility Check Passed"

    # Now optimize with all 1s objective
    @objective(separating_lp, Min, sum(x))
    optimize!(separating_lp)
    if is_solved_and_feasible(separating_lp)
        cut = get_cut_from_separating_solution(sepa, value.(x))
        push!(sepa.cutpool, cut)
    else
        @warn "All ones objective is unbounded"
    end

    # Now Optimize with p as objective
    @debug "Starting P* Optimization"
    p_star = argmin(point -> get_objective_value(point), translated_points)
    @objective(separating_lp, Min, sum(x[i] * p_star[i] for i in 1:dimension))
    optimize!(separating_lp)
    @debug "Finished P* Optimization"

    if is_solved_and_feasible(separating_lp)
        cut = get_cut_from_separating_solution(sepa, value.(x))
        push!(sepa.cutpool, cut)
    else
        error("p* as objective is unbounded")
    end

    if !is_EQ(scip, objective_value(separating_lp), 1.0)
        @error "Cannot find cut tight at p_star is not 1.0"
    end

    # Now transition to PRLP= 
    a_bar = value.(x)
    @constraint(separating_lp, sum(x[i] * p_star[i] for i in 1:dimension) == 1)
    set_time_limit_sec(separating_lp, 60.0)
    r_bar = filter(ray -> !is_zero(scip, dot(a_bar, ray)), rays)
    sort!(r_bar; by=ray -> abs(get_obj(get_generating_variable(ray))))
    objective_tried = 0
    marked = fill(false, length(r_bar))
    for ray in r_bar
        @objective(
            separating_lp, Min, sum(x[i] * ray[i] for i in 1:dimension)
        )
        optimize!(separating_lp)
        objective_tried += 1
        if primal_status(separating_lp) == MOI.FEASIBLE_POINT
            cut = get_cut_from_separating_solution(sepa, value.(x))
            push!(sepa.cutpool, cut)
        else
            if objective_tried > sepa.parameters.cut_limit
                break
            end
            println("Objective tried $(objective_tried)")
        end
        @debug "Generated $(length(sepa.cutpool)) cuts. Cutlimit is $(sepa.parameters.cut_limit)"
        if length(sepa.cutpool) >= sepa.parameters.cut_limit
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