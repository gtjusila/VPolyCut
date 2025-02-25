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
    time_limit = sepa.parameters.time_limit
    start_time = sepa.shared_data.start_time

    starting_lp_iter_count = SCIP.SCIPgetNLPIterations(scip)
    point_ray_collection = PointRayCollection()
    SCIP.SCIPstartProbing(scip)
    for (i, term) in enumerate(disjunction)
        @debug "Collecting Point and Ray from term $i"
        SCIP.SCIPnewProbingNode(scip)
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

            # Add point to point ray collection
            corner_point = project_point_to_nonbasic_space(
                nb_space, corner_polyhedron.lp_sol
            )
            add_point(
                point_ray_collection,
                corner_point
            )

            # Add rays to point ray collection
            # Rays have been complemented
            for ray in get_lp_rays(corner_polyhedron)
                add_ray(point_ray_collection, ray)
            end

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
