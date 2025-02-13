function gather_separating_solutions(
    prlp::PRLP,
    point_ray_collection::PointRayCollection;
    cut_limit::Int = typemax(Int),
    time_limit::Float64 = typemax(Float64),
    start_time::Float64 = time()
)
    separating_solutions = []
    # 6.1 The solution for feasibility is obtained from the callibration
    feasiblility = PRLPgetSolution(prlp)
    push!(separating_solutions, feasiblility)

    # 6.2 all ones
    @info "Trying all ones"
    all_ones = PRLPtryObjective(
        prlp, ones(SCIP.SCIP_Real, prlp.dimension); label = "all_ones"
    )
    if !isnothing(all_ones)
        push!(separating_solutions, all_ones)
    end

    # 6.3 Pstar 
    @info "Trying Pstar"
    pstar_index = argmin(1:length(get_points(point_ray_collection))) do i
        return get_objective_value(get_points(point_ray_collection)[i])
    end
    pstar = get_points(point_ray_collection)[pstar_index]
    pstar_separating_solution = PRLPtryObjective(
        prlp, as_dense_vector(pstar); label = "pstar"
    )
    if !isnothing(pstar_separating_solution)
        push!(separating_solutions, pstar_separating_solution)
    end

    # 6.4 Tighten Pstar and recalibrate PRLP
    @info "Tightening Pstar"
    PRLPtighten(prlp, pstar_index)
    @info "Recalibrating"
    if !PRLPcalibrate(prlp)
        throw(PStarNotTight())
    end
    push!(separating_solutions, PRLPgetSolution(prlp))

    # Now iterate over the rays 
    @info "Generating cuts from rays"
    a_bar = PRLPgetSolution(prlp)
    r_bar = filter(get_rays(point_ray_collection)) do ray
        return !is_zero(dot(a_bar, Vector(ray.coefficients))) # Use ray.coefficients instead of Ray to avoid type decution to AbstractArray
    end
    sort!(r_bar; by = ray -> abs(get_obj(get_generating_variable(ray))))
    objective_tried = 0

    for ray in r_bar
        r_run = PRLPtryObjective(
            prlp, as_dense_vector(ray);
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