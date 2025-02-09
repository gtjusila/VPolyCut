using SCIP
function get_point_ray_collection(
    scip::SCIP.SCIPData,
    disjunction::Disjunction
)
    point_ray_collection = PointRayCollection()

    @debug "Collecting Points and Rays from disjunction"
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
                throw("infeasible")
            end

            # Solve LP
            lp_feasible = solve_lp_relaxation(scip)
            if !lp_feasible
                @debug "Disjunctive term is detected infeasible during LP solve"
                throw("infeasible")
            end
            tableau = construct_tableau(scip)
            corner_polyhedron = construct_corner_polyhedron(tableau)

            # Add point to point ray collection
            add_point(
                point_ray_collection,
                get_lp_sol(corner_polyhedron),
                SCIP.SCIPgetLPObjval(scip),
                SCIP.SCIPgetSolOrigObj(scip, C_NULL)
            )

            # Add rays to point ray collection
            for ray in get_lp_rays(corner_polyhedron)
                add_ray(point_ray_collection, ray)
            end
        catch error
            if error != "infeasible"
                rethrow(error)
            end
            # if error is infeasible do nothing
        end
        SCIP.SCIPbacktrackProbing(scip, 0)
    end
    SCIP.SCIPendProbing(scip)
    return point_ray_collection
end