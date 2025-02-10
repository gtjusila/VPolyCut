function gather_separating_solutions(
    prlp::PRLP,
    point_ray_collection::PointRayCollection,
    non_basic_space::NonBasicSpace;
    cut_limit::Int = typemax(Int),
    time_limit::Float64 = typemax(Float64),
    start_time::Float64 = time(),
    scip::SCIP.SCIPData = SCIP.Optimizer().inner
)
    separating_solutions = []
    # 6.1 The solution for feasibility is obtained from the callibration
    feasiblility = PRLPgetSolution(prlp)
    push!(separating_solutions, feasiblility)

    # 6.2 all ones
    all_ones = PRLPtryObjective(
        prlp, ones(SCIP.SCIP_Real, prlp.dimension); label = "all_ones"
    )
    if !isnothing(all_ones)
        push!(separating_solutions, all_ones)
    end

    # 6.3 Pstar 
    pstar = argmin(
        x -> get_orig_objective_value(x), get_points(point_ray_collection)
    )
    pstar_projected = project_point_to_nonbasic_space(
        non_basic_space, get_point(pstar)
    )
    pstar_separating_solution = PRLPtryObjective(
        prlp, pstar_projected; label = "pstar"
    )
    if !isnothing(pstar_separating_solution)
        push!(separating_solutions, pstar_separating_solution)
    end

    # 6.4 Tighten Pstar and recalibrate PRLP
    p_star_zeroed = [!is_zero(scip, p) ? p : 0.0 for p in pstar_projected]
    PRLPtighten(prlp, p_star_zeroed)
    if !PRLPcalibrate(prlp)
        throw(PStarNotTight())
    end
    push!(separating_solutions, PRLPgetSolution(prlp))

    # Now iterate over the rays 
    a_bar = PRLPgetSolution(prlp)
    r_bar = filter(get_rays(point_ray_collection)) do ray
        projected_ray = project_ray_to_nonbasic_space(non_basic_space, ray)
        zeroed_ray = [!is_zero(scip, p) ? p : 0.0 for p in get_coefficients(projected_ray)]
        return !is_zero(scip, dot(a_bar, zeroed_ray))
    end
    sort!(r_bar; by = ray -> abs(get_obj(get_generating_variable(ray))))
    objective_tried = 0

    for ray in r_bar
        objective = get_coefficients(project_ray_to_nonbasic_space(non_basic_space, ray))
        r_run = PRLPtryObjective(
            prlp, objective;
            label = "prlp_$(objective_tried)"
        )
        objective_tried += 1
        if !isnothing(r_run)
            push!(separating_solutions, r_run)
            @debug "Generated $(length(separating_solutions)) cuts. Cutlimit is $(cut_limit)"
        end
        if objective_tried >= 2 * cut_limit
            break
        end
        if length(separating_solutions) >= cut_limit
            break
        end
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            break
        end
    end
    return separating_solutions
end