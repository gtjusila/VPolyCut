using SCIP

struct NodeInfeasible <: Exception end

"""
    get_point_ray_collection(scip::SCIP.SCIPData, disjunction::Disjunction)

Get the point ray collection for a given disjunction from the scip object.
Point ray collection will run in probing mode 
"""
function get_point_ray_collection(
    scip::SCIP.SCIPData,
    disjunction::Disjunction
)
    point_ray_collection = PointRayCollection()
    @debug "Collecting Points and Rays from disjunction"
    @info "$(memory_in_mb(point_ray_collection)) MB"
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
            tableau = construct_tableau(scip)
            corner_polyhedron = construct_corner_polyhedron(tableau)

            # Add point to point ray collection
            corner_point = get_solution_vector(tableau)
            add_point(
                point_ray_collection,
                corner_point,
                SCIP.SCIPgetLPObjval(scip),
                SCIP.SCIPgetSolOrigObj(scip, C_NULL)
            )

            # Add rays to point ray collection
            for ray in get_lp_rays(corner_polyhedron)
                add_ray(point_ray_collection, ray)
            end
        catch error
            if !(error isa NodeInfeasible)
                rethrow(error)
            end
            # if error is infeasible do nothing
        end
        SCIP.SCIPbacktrackProbing(scip, 0)
        @info "$(memory_in_mb(point_ray_collection)) MB"
    end
    SCIP.SCIPendProbing(scip)
    return point_ray_collection
end

function memory_in_mb(obj)
    bytes = Base.summarysize(obj)
    return bytes / (1024^2)  # Convert to MB
end