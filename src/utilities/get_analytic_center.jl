using SCIP
function get_analytic_center(scip::SCIP.SCIPData)
    SCIP.SCIPstartDive(scip)

    initial_algorithm = SCIP.get_parameter(scip, "lp/initalgorithm")
    resolve_algorithm = SCIP.get_parameter(scip, "lp/resolvealgorithm")
    check_dual_feasibility = SCIP.get_parameter(scip, "lp/checkdualfeas")
    disable_cutoff = SCIP.get_parameter(scip, "lp/disablecutoff")
    solution_polishing = SCIP.get_parameter(scip, "lp/solutionpolishing")
    check_stability = SCIP.get_parameter(scip, "lp/checkstability")
    check_farkas = SCIP.get_parameter(scip, "lp/checkfarkas")
    barrrier_convergence_tolerance = SCIP.get_parameter(scip, "numerics/barrierconvtol")

    SCIP.set_parameter(scip, "lp/initalgorithm", 'b')
    SCIP.set_parameter(scip, "lp/resolvealgorithm", 'b')
    SCIP.set_parameter(scip, "lp/checkdualfeas", 0)
    SCIP.set_parameter(scip, "lp/disablecutoff", 1)
    SCIP.set_parameter(scip, "lp/solutionpolishing", 0)
    SCIP.set_parameter(scip, "lp/checkstability", 0)
    SCIP.set_parameter(scip, "lp/checkfarkas", 0)
    SCIP.set_parameter(scip, "numerics/barrierconvtol", 1e-6)

    variable_pointers = collect_variable_pointers(scip)
    for var in variable_pointers
        SCIP.SCIPchgVarObjDive(scip, var, 2 * SCIP.SCIPvarGetObj(var))
    end
    lperror = Ref{SCIP.SCIP_Bool}(0)
    cutoff = Ref{SCIP.SCIP_Bool}(0)
    SCIP.SCIPsolveDiveLP(scip, -1, lperror, cutoff)
    if is_true(lperror[]) || is_true(cutoff[])
        throw(LPError())
    end
    sol = get_solution_vector(scip)
    SCIP.set_parameter(scip, "lp/initalgorithm", initial_algorithm)
    SCIP.set_parameter(scip, "lp/resolvealgorithm", resolve_algorithm)
    SCIP.set_parameter(scip, "lp/checkdualfeas", check_dual_feasibility)
    SCIP.set_parameter(scip, "lp/disablecutoff", disable_cutoff)
    SCIP.set_parameter(scip, "lp/solutionpolishing", solution_polishing)
    SCIP.set_parameter(scip, "lp/checkstability", check_stability)
    SCIP.set_parameter(scip, "lp/checkfarkas", check_farkas)
    SCIP.set_parameter(scip, "numerics/barrierconvtol", barrrier_convergence_tolerance)
    SCIP.SCIPendDive(scip)
    return sol
end