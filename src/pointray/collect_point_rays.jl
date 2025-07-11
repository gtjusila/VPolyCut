using SCIP

struct NodeInfeasible <: Exception end

"""
    get_point_ray_collection(scip::SCIP.SCIPData, disjunction::Disjunction)

Get the point ray collection for a given disjunction from the scip object.
Point ray collection will run in probing mode 
"""
function get_point_ray_collection(
    sepa::VPCSeparator
)
    scip = sepa.shared_data.scipd
    disjunction = sepa.shared_data.disjunction
    nb_space = sepa.shared_data.nonbasic_space
    time_limit = sepa.parameters.point_ray_collection_time_limit
    start_time = time()

    starting_lp_iter_count = SCIP.SCIPgetNLPIterations(scip)
    point_ray_collection = PointRayCollection()
    SCIP.SCIPstartProbing(scip)
    for (i, term) in enumerate(disjunction)
        @debug "Collecting Point and Ray from term $i"
        SCIP.SCIPnewProbingNode(scip)
        #Change objective
        @assert length(nb_space.variable_pointers) == length(nb_space.auxiliary_objective)
        for i in 1:length(nb_space.variable_pointers)
            SCIP.@SCIP_CALL SCIP.SCIPchgVarObjProbing(
                scip, nb_space.variable_pointers[i], nb_space.auxiliary_objective[i]
            )
        end
        try
            apply_bound_changes(term, scip)

            # Propegate
            prunable = propagate!(scip)
            if prunable
                @debug "Disjunctive term is detected infeasible during propagation"
                throw(NodeInfeasible())
            end

            # Solve LP
            lp_feasible = solve_lp_relaxation(scip)
            if !lp_feasible
                @debug "Disjunctive term is detected infeasible during LP solve"
                throw(NodeInfeasible())
            end
            #tableau = construct_tableau(scip)
            #corner_polyhedron = construct_corner_polyhedron(tableau)
            corner_polyhedron = CornerPolyhedron(scip, nb_space)
            push!(sepa.shared_data.original_points, corner_polyhedron.lp_sol)
            # Add point to point ray collection
            corner_point = project_point_to_nonbasic_space(
                nb_space, corner_polyhedron.lp_sol
            )
            add_point(
                point_ray_collection,
                corner_point
            )
            @debug "Nonzeros in corner point: $(SparseArrays.nnz(corner_point))"
            # Add rays to point ray collection
            # Rays have been complemented
            total_nonzeros = 0
            for ray in get_lp_rays(corner_polyhedron)
                total_nonzeros += SparseArrays.nnz(ray)
                add_ray(point_ray_collection, ray)
            end
            @debug "Total nonzeros in rays: $total_nonzeros"

            if time_limit < time() - start_time
                throw(TimeLimitExceededCollection())
            end
        catch error
            if !(error isa NodeInfeasible)
                rethrow(error)
            end
            # if error is infeasible do nothing
        end
        SCIP.SCIPbacktrackProbing(scip, 0)
    end
    SCIP.SCIPendProbing(scip)
    end_lp_iter_count = SCIP.SCIPgetNLPIterations(scip)
    sepa.statistics.point_ray_collection_lp_iterations =
        end_lp_iter_count - starting_lp_iter_count
    return point_ray_collection
end
