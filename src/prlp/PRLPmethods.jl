using SCIP
Base.show(io::IO, prlp::PRLP) = print(io, "PRLP with dimension $(prlp.dimension)")
Base.show(io::IO, scip::SCIP.SCIPData) = print(io, "SCIPData")

"""
    PRLPsetOutputStream(prlp::PRLP, stream::IO)

Set the output stream for the PRLP. By default, the output stream is set to `stdout`.
"""
function PRLPsetOutputStream(prlp::PRLP, stream::IO)
    prlp.output_stream = stream
end

"""
    PRLPinvalidateSolution(prlp::PRLP)

Invalidate the solution of the PRLP. This is useful if the objective is changed.
"""
function PRLPinvalidateSolution(prlp::PRLP)
    if prlp.solution_available
        prlp.solution_available = false
        prlp.solution_vector = []
    end
end

"""
    PRLPinvalidate(prlp::PRLP)

If PRLP is constructed, frees the LP and sets `lp_constructed` to false. Otherwise, do nothing.
"""
function PRLPinvalidate(prlp::PRLP)
    if prlp.lp_constructed
        if prlp.solution_available
            PRLPinvalidateSolution(prlp)
        end
        SCIP.SCIPlpiFree(address(prlp.lpi))
        prlp.lp_constructed = false
    end
end

"""
    PRLPaddPoint(prlp::PRLP, point::RealVector)

Add a single point to PRLP. If the LP for the PRLP have been constructed, this will invalidate the LP. Point will only be added if it is nonzero. 
"""
function PRLPaddPoint(prlp::PRLP, point::Point)
    @assert length(point) == prlp.dimension "Point dimension does not match PRLP dimension. Expected $(prlp.dimension), got $(length(point))"
    if (nnz(point) > 0)
        push!(prlp.points, point)
    else
        @warn "Trying to add zero point. Point is ignored."
    end
    PRLPinvalidate(prlp)
end

"""
    PRLPsetObjectiveLimit(prlp::PRLP, limit::SCIP.SCIP_Real)

Set the objective limit for the PRLP. This is useful if we want to stop the solve early if the objective value is greater than the limit.
"""
function PRLPsetObjectiveLimit(prlp::PRLP, limit::SCIP.SCIP_Real)
    SCIP.@SCIP_CALL SCIP.SCIPlpiSetRealpar(
        prlp.lpi,
        SCIP.SCIP_LPPAR_OBJLIM,
        limit
    )
end

"""
    PRLPaddRay(prlp::PRLP, ray::Ray)

Add a single ray to PRLP. If the LP for the PRLP have been constructed, this will invalidate the LP. Ray will only be added if it is nonzero.
"""
function PRLPaddRay(prlp::PRLP, ray::Ray)
    @assert length(ray) == prlp.dimension "Ray dimension does not match PRLP dimension"
    if (nnz(ray) > 0)
        push!(prlp.rays, ray)
    else
        @warn "Trying to add zero ray. Ray is ignored."
    end
    PRLPinvalidate(prlp)
end

"""
    PRLPconstructLP(prlp::PRLP)

Construct the LP for the PRLP, that is, creates the necessary variable and constraints. The objective will default to 0 and objective sense is minimize.
"""
function PRLPconstructLP(prlp::PRLP)
    # Create SCIP LPI object
    prlp.lp_constructed = true
    prlp.lpi = CPtr(SCIP.SCIP_LPI)
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(
        address(prlp.lpi), C_NULL, "PRLP", SCIP.SCIP_OBJSEN_MINIMIZE
    )

    # Add variables
    SCIP.SCIPlpiAddCols(
        prlp.lpi,
        prlp.dimension,
        zeros(SCIP.SCIP_Real, prlp.dimension),
        -SCIP.SCIPlpiInfinity(prlp.lpi) * ones(SCIP.SCIP_Real, prlp.dimension),
        SCIP.SCIPlpiInfinity(prlp.lpi) * ones(SCIP.SCIP_Real, prlp.dimension),
        ["x_$(i)" for i in 1:(prlp.dimension)],
        0,
        C_NULL,
        C_NULL,
        C_NULL
    )

    # Add point constraints
    beg::Vector{Cint} = []
    ind::Vector{Cint} = []
    val::Vector{SCIP.SCIP_Real} = []
    first = true
    for point in prlp.points
        push!(beg, length(ind))
        index, value = findnz(point)
        index .-= Cint(1)
        if first
            first = false
        end
        append!(ind, index)
        append!(val, value)
    end

    SCIP.SCIPlpiAddRows(
        prlp.lpi,
        length(beg),
        ones(SCIP.SCIP_Real, length(beg)),
        ones(SCIP.SCIP_Real, length(beg)) * SCIP.SCIPlpiInfinity(prlp.lpi),
        ["p_$(i)" for i in 1:length(beg)],
        length(ind),
        pointer(beg),
        pointer(ind),
        pointer(val)
    )

    # Clear the arrays 
    empty!(beg)
    empty!(ind)
    empty!(val)

    # Add ray constraints
    for ray in prlp.rays
        push!(beg, length(ind))
        index, value = findnz(ray)
        index .-= 1
        append!(ind, index)
        append!(val, value)
    end

    SCIP.SCIPlpiAddRows(
        prlp.lpi,
        length(beg),
        zeros(SCIP.SCIP_Real, length(beg)),
        ones(SCIP.SCIP_Real, length(beg)) * SCIP.SCIPlpiInfinity(prlp.lpi),
        ["r_$(i)" for i in 1:length(beg)],
        length(ind),
        pointer(beg),
        pointer(ind),
        pointer(val)
    )

    # LPI settings
    SCIP.@SCIP_CALL SCIP.SCIPlpiSetIntpar(
        prlp.lpi,
        SCIP.SCIP_LPPAR_PRESOLVING,
        1
    )
    #SCIP.@SCIP_CALL SCIP.SCIPlpiSetIntpar(
    #    prlp.lpi,
    #    SCIP.SCIP_LPPAR_LPINFO,
    #    1
    #)
end

"""
    PRLPsetObjective(prlp::PRLP, cost_vector::RealVector)

Set the cost vector for the PRLP. The cost vector must have the same dimension as the PRLP. The LP must be constructed before setting the objective.
"""
function PRLPsetObjective(prlp::PRLP, cost_vector::RealVector)
    @assert length(cost_vector) == prlp.dimension "Objective vector dimension does not match PRLP dimension"
    if !prlp.lp_constructed
        throw("LP must be constructed before setting objective")
    end
    PRLPinvalidateSolution(prlp)
    SCIP.@SCIP_CALL SCIP.SCIPlpiChgObj(
        prlp.lpi, prlp.dimension, [Cint(i) for i in 0:(prlp.dimension - 1)], cost_vector
    )
end

"""
    PRLPsetTimeLimit(prlp::PRLP, time_limit::SCIP.SCIP_Real)

Set the time_limit for the PRLP Solve in seconds. The LP must be constructed before setting the time limit.
"""
function PRLPsetTimeLimit(prlp::PRLP, time_limit::SCIP.SCIP_Real)
    if !prlp.lp_constructed
        throw("LP must be constructed before setting time limit")
    end
    SCIP.@SCIP_CALL SCIP.SCIPlpiSetRealpar(
        prlp.lpi, SCIP.SCIP_LPPAR_LPTILIM, time_limit
    )
end

"""
    PRLPsolveWithBarrier(prlp::PRLP)

Solve the PRLP using the barrier method. The LP must be constructed before solving.
"""
function PRLPsolveWithBarrier(prlp::PRLP)
    if !prlp.lp_constructed
        throw("LP must be constructed before solving")
    end
    if is_false(SCIP.SCIPlpiHasBarrierSolve())
        throw("LP Solver does not support barrier solve")
    end
    start_time = time()
    SCIP.@SCIP_CALL SCIP.SCIPlpiSolveBarrier(prlp.lpi, false)
    end_time = time()
    prlp.last_solve_time = end_time - start_time
    prlp.last_simplex_iterations = -1

    # For barrier method we only get solution if LP is optimal
    if is_true(SCIP.SCIPlpiIsOptimal(prlp.lpi))
        prlp.solution_available = true
        prlp.solution_vector = zeros(SCIP.SCIP_Real, prlp.dimension)
        obj_val = Ref{SCIP.SCIP_Real}(0)
        SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(prlp.lpi, obj_val, pointer(prlp.solution_vector),
            C_NULL,
            C_NULL, C_NULL
        )
        prlp.solution_objective = obj_val[]
    end
end

"""
    PRLPsolveWithPrimalSimplex(prlp::PRLP)

Solve the PRLP using the primal simplex method. The LP must be constructed before solving.
"""
function PRLPsolveWithPrimalSimplex(prlp::PRLP)
    if !prlp.lp_constructed
        throw("LP must be constructed before solving")
    end
    if is_false(SCIP.SCIPlpiHasPrimalSolve())
        throw("LP Solver does not support primal simplex solve")
    end

    start_time = time()
    SCIP.@SCIP_CALL SCIP.SCIPlpiSolvePrimal(prlp.lpi)
    end_time = time()
    prlp.last_solve_time = end_time - start_time

    simplex_iteration = Ref{Cint}(0)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetIterations(prlp.lpi, simplex_iteration)
    prlp.last_simplex_iterations = simplex_iteration[]

    if is_true(SCIP.SCIPlpiIsPrimalFeasible(prlp.lpi))
        prlp.solution_available = true
        prlp.solution_vector = zeros(SCIP.SCIP_Real, prlp.dimension)
        obj_val = Ref{SCIP.SCIP_Real}(0)
        SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(prlp.lpi, obj_val, pointer(prlp.solution_vector),
            C_NULL,
            C_NULL, C_NULL
        )
        prlp.solution_objective = obj_val[]
    end
end

"""
    PRLPsolveWithDualSimplex(prlp::PRLP)

Solve the PRLP using the primal simplex method. The LP must be constructed before solving.
"""
function PRLPsolveWithDualSimplex(prlp::PRLP)
    if !prlp.lp_constructed
        throw("LP must be constructed before solving")
    end
    if is_false(SCIP.SCIPlpiHasDualSolve())
        throw("LP Solver does not support dual simplex solve")
    end

    start_time = time()
    SCIP.@SCIP_CALL SCIP.SCIPlpiSolveDual(prlp.lpi)
    end_time = time()
    prlp.last_solve_time = end_time - start_time

    simplex_iteration = Ref{Cint}(0)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetIterations(prlp.lpi, simplex_iteration)
    prlp.last_simplex_iterations = simplex_iteration[]

    if is_true(SCIP.SCIPlpiIsPrimalFeasible(prlp.lpi))
        prlp.solution_available = true
        prlp.solution_vector = zeros(SCIP.SCIP_Real, prlp.dimension)
        obj_val = Ref{SCIP.SCIP_Real}(0)
        SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(prlp.lpi, obj_val, pointer(prlp.solution_vector),
            C_NULL,
            C_NULL, C_NULL
        )
        prlp.solution_objective = obj_val[]
    end
end

"""
    PRLPsetSolvingAlgorithm(prlp::PRLP, algorithm::PRLPsolveAlgorithm)

Set the solving algorithm for the PRLP. 
"""
function PRLPsetSolvingAlgorithm(prlp::PRLP, algorithm::PRLPsolveAlgorithm)
    prlp.solve_algorithm = algorithm
end

"""
    PRLPsolve(prlp::PRLP)

Solve the PRLP using the specified algorithm. The LP must be constructed before solving.
"""
function PRLPsolve(prlp::PRLP)
    if prlp.solve_algorithm == PRIMAL_SIMPLEX
        PRLPsolveWithPrimalSimplex(prlp)
    elseif prlp.solve_algorithm == BARRIER
        PRLPsolveWithBarrier(prlp)
    elseif prlp.solve_algorithm == DUAL_SIMPLEX
        PRLPsolveWithDualSimplex(prlp)
    else
        throw("Unknown solving algorithm")
    end
end
"""
    PRLPisSolutionAvailable(prlp::PRLP)

Check if the solution is available for the PRLP.
"""
function PRLPisSolutionAvailable(prlp::PRLP)
    return prlp.solution_available
end

"""
    PRLPgetSolution(prlp::PRLP)

Get the solution of the PRLP. Make sure solution is available before calling this function, otherwise this will throw an error.
"""
function PRLPgetSolution(prlp::PRLP)
    if !prlp.solution_available
        throw("Solution is not available")
    end
    return prlp.solution_vector
end

"""
    PRLPgetObjective(prlp::PRLP)

Get the objective value of the solution of the PRLP. Make sure solution is available before calling this function, otherwise this will throw an error.
"""
function PRLPgetObjective(prlp::PRLP)
    if !prlp.solution_available
        throw("Solution is not available")
    end
    return prlp.solution_objective
end

"""
    PRLPtighten(prlp::PRLP, point::RealVector)

Set the inequality constraint corresponding to point to be an equality constraint.
"""
function PRLPtighten(prlp::PRLP, index::Int)
    if prlp.lp_constructed == false
        throw("LP must be constructed before tightening")
    end
    SCIP.SCIPlpiChgSides(prlp.lpi, 1, [Cint(index - 1)], [1.0], [1.0])
end

"""
    PRLPgetLastSolveTime(prlp::PRLP)

Get the time taken to solve the last PRLP.
"""
function PRLPgetLastSolveTime(prlp::PRLP)
    return prlp.last_solve_time
end

"""
    PRLPgetLastSimplexIterations(prlp::PRLP)

Get the number of simplex iterations taken to solve the last PRLP. This only give meaningful result if the last solve was done with either primal or dual simplex.
For barrier solve it will be set to -1.
"""
function PRLPgetLastSimplexIterations(prlp::PRLP)
    return prlp.last_simplex_iterations
end

"""
    PRLPdestroy(prlp::PRLP) 

Destroy the PRLP and free the memory.
"""
function PRLPdestroy(prlp::PRLP)
    PRLPinvalidate(prlp)
end

"""
    PRLPgetSolvingAlgorithm(prlp::PRLP)

Get the solving algorithm of the PRLP.
"""
function PRLPgetSolvingAlgorithm(prlp::PRLP)
    return prlp.solve_algorithm
end

"""
    PRLPgetLastTerminationStatus(prlp::PRLP)

Get the termination status of the last solve of the PRLP.
"""
function PRLPgetLastTerminationStatus(prlp::PRLP)
    return lpi_termination_status(prlp.lpi)
end

function lpi_termination_status(lpi::CPtr{SCIP.SCIP_LPI})::TerminationStatus
    if is_false(SCIP.SCIPlpiWasSolved(lpi))
        return LPI_NOT_SOLVED
    end
    if is_true(SCIP.SCIPlpiIsOptimal(lpi))
        return LPI_OPTIMAL
    end
    if is_true(SCIP.SCIPlpiIsTimelimExc(lpi))
        return LPI_TIME_LIMIT_EXCEEDED
    end
    if is_true(SCIP.SCIPlpiExistsPrimalRay(lpi))
        return LPI_UNBOUNDED
    end
    @error "Undandled termination status. Status code $(SCIP.SCIPlpiGetInternalStatus(lpi))"
    return LPI_NOT_SOLVED
end

"""
    PRLPtryObjective(prlp::PRLP, objective::RealVector; time_limit = 30.0, label::String = "")

Try to optimize the given objective. If the solution is found, return the solution vector, otherwise return nothing.
"""
function PRLPtryObjective(
    prlp::PRLP, objective::ObjectiveFunction
)::Union{Vector{SCIP.SCIP_Real},Nothing}
    # Warm start if 1. last solve is not good (otherwise LP solver could warmstart on its own)
    # 2. Objective is savable and 3. Warm start is allowed
    if !prlp.last_solve_good && objective.savable && prlp.allow_warm_start
        @debug "Last solve was not good. Trying to restore basis"
        PRLPrestoreBasis(prlp)
    end

    # Set time limit, objective and try the solve
    PRLPsetTimeLimit(prlp, objective.time_limit)
    PRLPsetObjective(prlp, objective.direction)
    prlp.last_solve_good = false
    PRLPsolve(prlp)
    PRLPrecordSolveStats(prlp, objective.label)

    # If lpi is optimal then solve is good, save the basis
    if is_true(SCIP.SCIPlpiIsOptimal(prlp.lpi))
        prlp.last_solve_good = true
        # if warmstart is allowed and basis is savable, save the basis
        # test only store first basis
        if objective.savable && prlp.allow_warm_start &&
            objective.label == "pstar_feasibility"
            PRLPstoreBasis(prlp)
        end
    end

    # Return the solution if available
    if PRLPisSolutionAvailable(prlp)
        return PRLPgetSolution(prlp)
    end

    # If solution is not available, check if on_fail is set
    if !isnothing(objective.on_fail)
        throw(objective.on_fail)
    end

    return nothing
end

"""
    PRLPrecordSolveStats(prlp::PRLP, label::String)

Record the solve statistics of the PRLP. This will record the solve algorithm, solve time, simplex iterations, primal status and termination status.
"""
function PRLPrecordSolveStats(prlp::PRLP, label::String)
    stats = Dict(
        "solve_algorithm" => PRLPgetSolvingAlgorithm(prlp),
        "solve_time" => PRLPgetLastSolveTime(prlp),
        "simplex_iterations" => PRLPgetLastSimplexIterations(prlp),
        "primal_status" => PRLPisSolutionAvailable(prlp),
        "termination_status" => PRLPgetLastTerminationStatus(prlp),
        "objective_name" => label
    )
    push!(prlp.solve_statistics, stats)
end

"""
    PRLPstoreBasis(prlp::PRLP)

Store the basis of the PRLP. This is useful if we want to warm start the LP solver.
"""
function PRLPstoreBasis(prlp::PRLP)
    blkmem = SCIP.SCIPblkmem(prlp.scip)
    PRLPfreeBasis(prlp)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetState(prlp.lpi, blkmem, prlp.last_good_state)
    @assert !(prlp.last_good_state[] == C_NULL) "Failed to store basis"
end

"""
    PRLPrestoreBasis(prlp::PRLP)

Restore the last known optimal basis of the PRLP.
"""
function PRLPrestoreBasis(prlp::PRLP)
    if prlp.last_good_state[] != C_NULL
        blkmem = SCIP.SCIPblkmem(prlp.scip)
        prlp.n_basis_restart += 1
        SCIP.@SCIP_CALL SCIP.SCIPlpiSetState(prlp.lpi, blkmem, prlp.last_good_state[])
    end
end

"""
    PRLPfreeBasis(prlp::PRLP)

Free the basis of the PRLP. This is useful if we want to free the memory.
"""
function PRLPfreeBasis(prlp::PRLP)
    blkmem = SCIP.SCIPblkmem(prlp.scip)
    if prlp.last_good_state[] != C_NULL
        SCIP.SCIPlpiFreeState(prlp.lpi, blkmem, prlp.last_good_state)
    end
end