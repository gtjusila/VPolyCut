function solve_separation_subproblems(sepa::VPCSeparator)
    scip = sepa.scipd

    lp_solution = get_solution_vector(sepa.complemented_tableau)
    lp_solution = project(sepa.projection, lp_solution)
    dimension = length(lp_solution)

    separating_lp = Model(HiGHS.Optimizer)
    set_optimizer_attribute(separating_lp, "output_flag", false)
    zeroed = 0
    # translate points to non-basic space
    points = get_points(sepa.point_ray_collection)
    translated_points = map(points) do point
        return CornerPoint(get_point(point) - lp_solution, get_objective_value(point))
    end

    rays = get_rays(sepa.point_ray_collection)

    vars = map(1:get_nvars(sepa.complemented_tableau)) do i
        return get_var_from_column(sepa.complemented_tableau, i)
    end
    c_bar = map(vars) do var
        return get_obj(var)
    end
    lp_sol = map(vars) do var
        return get_sol(var)
    end
    lp_obj = dot(c_bar, lp_sol)
    @debug "LP Objective is $(lp_obj) SepaLP OBJ is $(sepa.lp_obj)"
    # Construct PRLP0
    @variable(separating_lp, x[1:dimension] >= 0)
    for point in translated_points
        vertex = get_point(point)
        vertex = map(vertex) do x
            if is_zero(sepa.scipd, x)
                if x != 0
                    zeroed += 1
                end
                return 0.0
            else
                return x
            end
        end
        @constraint(separating_lp, sum(x[i] * vertex[i] for i in 1:dimension) >= 1)
    end

    for ray in rays
        coefficient = get_coefficients(ray)
        coefficient = map(coefficient) do x
            if is_zero(sepa.scipd, x)
                if x != 0
                    zeroed += 1
                end
                return 0.0
            else
                return x
            end
        end
        @constraint(separating_lp, sum(x[i] * coefficient[i] for i in 1:dimension) >= 0)
    end

    @debug "Zeroed $(zeroed) entries"

    # Check if disjunctive lower bound is worse than optimal LP objective
    disjunctive_lower_bount = minimum(x -> get_objective_value(x), translated_points)
    @debug "Disjunctive Lower Bound is: $(disjunctive_lower_bount) and LP Objective is: $(sepa.lp_obj)"
    if !is_GT(scip, disjunctive_lower_bount, sepa.lp_obj)
        throw(FailedDisjunctiveLowerBoundTest())
    end

    # check if the LP is feasible by optimizing using all 0s objective
    @objective(separating_lp, Min, 0)
    write_to_file(separating_lp, "separating_lp.lp")
    println("Starting Check Feasibility")
    set_time_limit_sec(separating_lp, 120.0)
    optimize!(separating_lp)
    if primal_status(separating_lp) != MOI.FEASIBLE_POINT
        throw(FailedToProvePRLPFeasibility())
    end
    println("Feasibility Check Passed")
    # Now optimize with all 1s objective
    @objective(separating_lp, Min, sum(x))
    optimize!(separating_lp)
    if is_solved_and_feasible(separating_lp)
        cut = get_cut_from_separating_solution(sepa, value.(x))
        push!(sepa.cutpool, cut)
    else
        error("All ones objective is unbounded")
    end

    # Now Optimize with p as objective
    println("Starting P* Optimization")
    p_star = argmin(point -> get_objective_value(point), translated_points)
    @objective(separating_lp, Min, sum(x[i] * p_star[i] for i in 1:dimension))
    optimize!(separating_lp)
    println("Finished P* Optimization")

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
    set_time_limit_sec(separating_lp, 10.0)
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
        if is_solved_and_feasible(separating_lp)
            cut = get_cut_from_separating_solution(sepa, value.(x))
            push!(sepa.cutpool, cut)
        else
            if objective_tried > sepa.parameters.cut_limit * 2
                break
            end
            println("Objective tried $(objective_tried)")
        end
        println(
            "Generated $(length(sepa.cutpool)) cuts. Cutlimit is $(sepa.parameters.cut_limit)"
        )
        if length(sepa.cutpool) >= sepa.parameters.cut_limit
            break
        end
    end
    add_all_cuts!(sepa.cutpool, sepa)
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