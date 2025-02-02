using Xpress
using SCIP

function solve_separation_subproblems_with_scip(sepa::VPCSeparator)
    # Setup aliases
    scip = sepa.scipd

    # Get current LP Solution projected onto the non-basic variables.
    # This will be the origin point 
    current_lp_solution = get_solution_vector(sepa.complemented_tableau)
    projected_current_lp_solution = project(sepa.projection, current_lp_solution)
    problem_dimension = length(projected_current_lp_solution)

    # First, translate the points so that the lp_solution is at the origin 
    projected_point_collection = get_points(sepa.point_ray_collection)
    point_collection_in_nonbasic_space = map(projected_point_collection) do corner_point
        return CornerPoint(
            get_point(corner_point) - projected_current_lp_solution,
            get_objective_value(corner_point),
            get_orig_objective_value(corner_point)
        )
    end
    ray_collection = get_rays(sepa.point_ray_collection)

    # Create SCIP problem
    subscip = CPtr(SCIP.SCIP_)
    SCIP.@SCIP_CALL SCIP.SCIPcreate(address(subscip))
    SCIP.@SCIP_CALL SCIP.SCIPincludeDefaultPlugins(subscip)
    SCIP.@SCIP_CALL SCIP.SCIPcreateProbBasic(subscip, "PRLP")
    SCIP.@SCIP_CALL SCIP.SCIPsetSubscipsOff(subscip, 0)
    SCIP.@SCIP_CALL SCIP.SCIPenableReoptimization(subscip, true)
    # Variables 
    vars = [CPtr(SCIP.SCIP_Var) for i in 1:problem_dimension]
    for i in 1:problem_dimension
        SCIP.@SCIP_CALL SCIP.SCIPcreateVarBasic(
            subscip, address(vars[i]), "x_$(i)",
            -SCIP.SCIPinfinity(subscip), SCIP.SCIPinfinity(subscip),
            0, SCIP.SCIP_VARTYPE_CONTINUOUS
        )
        SCIP.@SCIP_CALL SCIP.SCIPaddVar(subscip, vars[i])
    end
    # Point constraint
    n_points = length(point_collection_in_nonbasic_space)
    point_constraints = [CPtr(SCIP.SCIP_Cons) for i in 1:n_points]

    for i in 1:n_points
        vertex = get_point(point_collection_in_nonbasic_space[i])
        SCIP.@SCIP_CALL SCIP.SCIPcreateConsBasicLinear(
            subscip, address(point_constraints[i]), "p_$(i)",
            0, C_NULL, C_NULL, 1, SCIP.SCIPinfinity(scip)
        )

        for j in 1:problem_dimension
            SCIP.@SCIP_CALL SCIP.SCIPaddCoefLinear(
                subscip, point_constraints[i], vars[j], vertex[j]
            )
        end
        SCIP.@SCIP_CALL SCIP.SCIPaddCons(subscip, point_constraints[i])
    end
    # Ray constraint
    n_rays = length(ray_collection)
    ray_constraints = [CPtr(SCIP.SCIP_Cons) for i in 1:n_rays]
    for i in 1:n_rays
        ray = get_coefficients(ray_collection[i])
        SCIP.@SCIP_CALL SCIP.SCIPcreateConsBasicLinear(
            subscip, address(ray_constraints[i]), "r_$(i)",
            0, C_NULL, C_NULL, 0, SCIP.SCIPinfinity(scip)
        )

        for j in 1:problem_dimension
            SCIP.@SCIP_CALL SCIP.SCIPaddCoefLinear(
                subscip, ray_constraints[i], vars[j], ray[j]
            )
        end

        SCIP.@SCIP_CALL SCIP.SCIPaddCons(subscip, ray_constraints[i])
    end

    # Use P* as objective
    p_star = argmin(point -> get_objective_value(point), point_collection_in_nonbasic_space)
    for i in 1:problem_dimension
        SCIP.@SCIP_CALL SCIP.SCIPchgVarObj(subscip, vars[i], p_star[i])
    end
    SCIP.@SCIP_CALL SCIP.SCIPwriteOrigProblem(subscip, "separating_test.lp", C_NULL, true)
    SCIP.@SCIP_CALL SCIP.SCIPsolve(subscip)

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

function get_solve_stat(model::JuMP.AbstractModel, obj_name::String)
    return Dict(
        "solve_time" => solve_time(model),
        #"simplex_iterations" => simplex_iterations(model),
        "primal_status" => primal_status(model),
        "dual_status" => dual_status(model),
        "termination_status" => termination_status(model),
        "objective_name" => obj_name
    )
end
