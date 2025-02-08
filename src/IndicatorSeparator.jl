using SCIP

"""
Indicator Separator
# Fields
- scipd::SCIP.SCIPData Reference to the SCIPData object
"""
@kwdef mutable struct IndicatorSeparator <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
    called::Int = 0
end

function SCIP.exec_lp(sepa::IndicatorSeparator)
    ## Warning Spaghetti Code Ahead ##

    # Hide stdout 
    #temp = stdout
    #redirect_stdout(devnull)

    # Aliasing for easier call
    scip = sepa.scipd

    # Get the indicator variable dictionary
    indicator_var_dict = getIndicatorVarDict(scip)

    # A dict to store variable fixings found during the process
    fixing = Dict{Ptr{SCIP.SCIP_VAR},SCIP.SCIP_Real}()
    # Only do first 100
    cnt = 0

    # Loop through the indicator variables
    for key in keys(indicator_var_dict)
        entry = indicator_var_dict[key]

        # For now we only process indicator variables that are fractional or not satisfied
        if is_indicator_satisfied(scip, entry)
            continue
        end

        # Some prints incase we need to debug in the future
        println("Var: ", unsafe_string(SCIP.SCIPvarGetName(key)))
        println("Is Negated: ", entry.is_negated)
        println("Current Value: ", SCIP.SCIPvarGetLPSol(key))
        #=
        println(
            "Negated Var: ",
            unsafe_string(
                SCIP.SCIPvarGetName(entry.negated_var)
            )
        )
        =#
        println(
            "Constraints: ",
            [
                unsafe_string(SCIP.SCIPconsGetName(cons)) for
                cons in entry.constraints
            ]
        )
        println(
            "Negated Constraints: ",
            [
                unsafe_string(SCIP.SCIPconsGetName(cons)) for
                cons in entry.negated_constraints
            ]
        )
        println("Is satisfied? ", is_indicator_satisfied(scip, entry))
        # End of prints
        #=
                scip_tableau = construct_tableau_with_constraint_matrix(sepa.scipd)
                complemented_tableau = ComplementedTableau(scip_tableau)
                projection = create_projection_to_nonbasic_space(complemented_tableau)
                point_ray_collection = PointRayCollection(scip; projection=projection)
                sol1 = CPtr(SCIP.SCIP_Sol)
                SCIP.@SCIP_CALL SCIP.SCIPcreateLPSol(scip, address(sol1), C_NULL)
              =#
        # Start Probing model
        SCIP.@SCIP_CALL SCIP.SCIPstartProbing(scip)
        SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(scip)

        # Fix Binary to 1
        SCIP.@SCIP_CALL SCIP.SCIPfixVarProbing(scip, key, 1.0)
        for cons in entry.constraints
            slack = SCIP.SCIPgetSlackVarIndicator(cons)
            SCIP.@SCIP_CALL SCIP.SCIPfixVarProbing(scip, slack, 0.0)
        end
        for cons in entry.negated_constraints
            slack = SCIP.SCIPgetSlackVarIndicator(cons)
            SCIP.@SCIP_CALL SCIP.SCIPchgVarUbProbing(scip, slack, SCIP.SCIPinfinity(scip))
        end

        if propagate!(scip)
            # Prunable
            println("Prunable")
            fixing[key] = 0.0
            SCIP.SCIPendProbing(scip)
            continue
        end

        if !solve_lp_relaxation(scip)
            # Infeasible LP
            println("Infeasible LP")
            fixing[key] = 0.0
            SCIP.SCIPendProbing(scip)
            continue
        end
        # Done Now recover corner CornerPolyhedron
        #=
        # Get Optimal Tableau
        left_tableau = construct_tableau(scip)

        # Complement the same columns as the root tableau
        complemented_columns = get_complemented_columns(complemented_tableau)
        for i in complemented_columns
            var = get_var_from_column(left_tableau, i)
            complement_column(left_tableau, var)
        end

        left_corner_polyhedron = construct_corner_polyhedron(left_tableau)
        left_objective_value = SCIP.SCIPgetLPObjval(scip)
        left_orig_objective_value = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
        left_basic_solution = get_lp_sol(left_corner_polyhedron)

        add_point(
            point_ray_collection, left_basic_solution, left_objective_value,
            left_orig_objective_value
        )
        for ray in get_lp_rays(left_corner_polyhedron)
            add_ray(point_ray_collection, ray)
        end
        sol2 = CPtr(SCIP.SCIP_Sol)
        SCIP.@SCIP_CALL SCIP.SCIPcreateLPSol(scip, address(sol2), C_NULL)
        =#
        #Done with the first corner polyhedron
        # Again now some prints for debug
        ##println("Objective After Fixing Binary to 1: $(SCIP.SCIPgetLPObjval(scip))")
        #println("Constraint is now satisfied? ", is_indicator_satisfied(scip, entry))
        @assert is_indicator_satisfied(scip, entry)
        SCIP.SCIPbacktrackProbing(scip, 0)

        # Now fixing binary to 0 we call this right fixing
        SCIP.SCIPnewProbingNode(scip)
        SCIP.@SCIP_CALL SCIP.SCIPfixVarProbing(scip, key, 0.0)
        for cons in entry.negated_constraints
            slack = SCIP.SCIPgetSlackVarIndicator(cons)
            SCIP.@SCIP_CALL SCIP.SCIPfixVarProbing(scip, slack, 0.0)
        end
        for cons in entry.constraints
            slack = SCIP.SCIPgetSlackVarIndicator(cons)
            SCIP.@SCIP_CALL SCIP.SCIPchgVarUbProbing(scip, slack, SCIP.SCIPinfinity(scip))
        end

        if propagate!(scip)
            # Prunable
            fixing[key] = 1.0
            println("Prunable")
            SCIP.SCIPendProbing(scip)
            continue
        end

        if !solve_lp_relaxation(scip)
            # Infeasible LP
            fixing[key] = 1.0
            println("Infeasible LP")
            SCIP.SCIPendProbing(scip)
            continue
        end
        SCIP.@SCIP_CALL SCIP.SCIPendProbing(scip)
        # Done Now recover corner CornerPolyhedron
        #=
         # Get Optimal Tableau
         right_tableau = construct_tableau(scip)

         # Complement the same columns as the root tableau
         complemented_columns = get_complemented_columns(complemented_tableau)
         for i in complemented_columns
             var = get_var_from_column(right_tableau, i)
             complement_column(right_tableau, var)
         end

         right_corner_polyhedron = construct_corner_polyhedron(right_tableau)
         right_objective_value = SCIP.SCIPgetLPObjval(scip)
         right_orig_objective_value = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
         right_basic_solution = get_lp_sol(right_corner_polyhedron)

         add_point(
             point_ray_collection, right_basic_solution, right_objective_value,
             right_orig_objective_value
         )
         for ray in get_lp_rays(right_corner_polyhedron)
             add_ray(point_ray_collection, ray)
         end
         sol3 = CPtr(SCIP.SCIP_Sol)
         SCIP.@SCIP_CALL SCIP.SCIPcreateLPSol(scip, address(sol3), C_NULL)
         # Done with the second corner polyhedron

         println("Objective After Fixing Binary to 0: $(SCIP.SCIPgetLPObjval(scip))")
         @assert is_indicator_satisfied(scip, entry)
         SCIP.SCIPendProbing(scip)

         # Now we generate cut
         points = get_points(point_ray_collection)
         rays = get_rays(point_ray_collection)
         original_lp_solution = get_solution_vector(complemented_tableau)
         projected_original_lp_solution = project(projection, original_lp_solution)

         points_in_nonbasic_space = [
             get_point(p) - projected_original_lp_solution for p in points
         ]
         rays_in_nonbasic_space = [
             get_coefficients(r) for r in rays
         ]

         # Construct the PRLP
         prlp = PRLP(length(projected_original_lp_solution))
         for point in points_in_nonbasic_space
             p = [!is_zero(scip, x) ? x : 0.0 for x in point]
             PRLPaddPoint(prlp, p)
         end
         for ray in rays_in_nonbasic_space
             r = [!is_zero(scip, x) ? x : 0.0 for x in ray]
             PRLPaddRay(prlp, r)
         end

         function convert_cut_from_separating_vector(vector)
             lp_solution = get_solution_vector(complemented_tableau)
             lp_solution = get_uncomplemented_vector(lp_solution, complemented_tableau)

             separating_solution = undo_projection(projection, vector)
             separating_solution = get_uncomplemented_vector(
                 separating_solution, complemented_tableau
             )
             b = dot(separating_solution, lp_solution) + 1

             cut_vector, b = convert_standard_inequality_to_general(
                 scip, complemented_tableau, separating_solution, b
             )

             row = add_sepa_row!(
                 scip,
                 sepa,
                 -cut_vector,
                 get_problem_variables_pointers(complemented_tableau),
                 -b
             )
             println("Hello")
             SCIP.@SCIP_CALL SCIP.SCIPprintRow(scip, row[], C_NULL)
             SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip, row)
         end

         PRLPconstructLP(prlp)
         PRLPsetSolvingAlgorithm(prlp, PRIMAL_SIMPLEX)
         PRLPsetTimeLimit(prlp, 60.0)

         #for point in points_in_nonbasic_space
         #    p = [!is_zero(scip, x) ? x : 0.0 for x in point]
         #    PRLPsetObjective(prlp, p)
         #    PRLPsolve(prlp)
         #    if PRLPisSolutionAvailable(prlp)
         #        sol = PRLPgetSolution(prlp)
         #        convert_cut_from_separating_vector(sol)
         #    end
         #end

         for point in points_in_nonbasic_space
             p = [!is_zero(scip, x) ? x : 0.0 for x in point]
             PRLPtighten(prlp, p)
         end
         PRLPsetObjective(
             prlp, zeros(SCIP.SCIP_Real, length(projected_original_lp_solution))
         )
         PRLPsolve(prlp)
         if PRLPisSolutionAvailable(prlp)
             sol = PRLPgetSolution(prlp)
             convert_cut_from_separating_vector(sol)
         end
         cnt += 1
         =#
    end

    @info "Number of indicator $(length(keys(indicator_var_dict)))"
    fixed_cnt = 0
    for key in keys(fixing)
        if is_EQ(scip, fixing[key], 1.0) && is_EQ(scip, SCIP.SCIPvarGetLbLocal(key), 1.0)
            continue
        end
        if is_EQ(scip, fixing[key], 0.0) && is_EQ(scip, SCIP.SCIPvarGetUbLocal(key), 0.0)
            continue
        end
        infeasible = Ref{SCIP.SCIP_Bool}()
        fixed = Ref{SCIP.SCIP_Bool}()
        SCIP.@SCIP_CALL SCIP.SCIPfixVar(
            scip, key, fixing[key], infeasible, fixed
        )
        if is_true(infeasible[])
            @warn "Infeasible fixing"
        end
        if is_false(fixed[])
            @error "Failed to fix"
        end
        fixed_cnt = 1
    end
    #redirect_stdout(temp)
    @info "Number of fixing $(fixed_cnt)"
    if fixed_cnt > 0
        return SCIP.SCIP_REDUCEDDOM
    end
    return SCIP.SCIP_DIDNOTFIND
end

function include_indicator_sepa(scip::SCIP.SCIPData)
    sepa = IndicatorSeparator(; scipd=scip)
    SCIP.include_sepa(scip.scip[], scip.sepas, sepa; freq=0, usessubscip=true)
end