#
# This file contains the exec_lp function for VPC separation
#
using SCIP
using JuMP
using LinearAlgebra
import MathOptInterface as MOI

"""
VPolyhedral Cut Separator

Implementation of Algorithm 1 from 
Balas, Egon, and Aleksandr M. Kazachkov. "V-polyhedral disjunctive cuts." arXiv preprint arXiv:2207.13619 (2022).
Disjunction are obtained from partial branch and bound trees
"""

# The first layer is a wrapper for logging
function SCIP.exec_lp(sepa::VPCSeparator)
    # Wrapper for writing logs
    file_logger = setup_file_logger(
        joinpath(sepa.parameters.log_directory, "vpc_separation.log")
    )
    console_logger = ConsoleLogger(Logging.Debug)
    logger = TeeLogger(file_logger, console_logger)

    with_logger(logger) do
        return _exec_lp(sepa)
    end
end

# The second layer is a try block to catch potential error arising during cut generation
function _exec_lp(sepa::VPCSeparator)
    # Aliasing for easier call
    scip = sepa.scipd

    # Initialization: start timer, set separated to false, and increment call count 
    @info "VPC Separator Called"
    sepa.statistics.called += 1

    @info "Root node LP took $(SCIP.SCIPgetNRootFirstLPIterations(scip)) iterations"
    # Check Preconditions and handle accordingly
    if sepa.should_be_skipped
        return SCIP.SCIP_DIDNOTRUN
    end
    if !is_numeric_scip_set()
        set_numeric_scip(scip)
    end
    if SCIP.SCIPgetStage(scip) != SCIP.SCIP_STAGE_SOLVING
        return SCIP.SCIP_DIDNOTRUN
    end
    if SCIP.SCIPisLPSolBasic(scip) == 0
        # LP Solution is not basic 
        return SCIP.SCIP_DELAYED
    end
    if SCIP.SCIPgetLPSolstat(scip) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        # LP Solution is not optimal 
        return SCIP.SCIP_DELAYED
    end
    @info "Finished checking necessary Preconditions"

    # Do everything in a try block to ensure time limit requirement
    @info "Starting separation subroutine"
    try
        # Call Separation Subroutine
        vpolyhedralcut_separation(sepa)
    catch error
        # If error occured while in probing mode, we need to cleanup
        if is_true(SCIP.SCIPinProbing(scip))
            SCIP.SCIPendProbing(scip)
        end

        # Mark sepa to be skipped on next call
        # Any errors happening means sepa is not good for the current problem
        sepa.should_be_skipped = true

        # Handle different types of errors
        if error isa TimeLimitExceededCollection
            @debug "Time Limit Exceeded Collection"
            sepa.termination_status = TIME_LIMIT_EXCEEDED_COLLECTION
        elseif error isa LPError
            @debug "LP error"
            sepa.termination_status = LP_ERROR
        elseif error isa TimeLimitExceededBranchAndBound
            @debug "Time Limit Exceeded Branch and Bound"
            sepa.termination_status = TIME_LIMIT_EXCEEDED_BRANCHANDBOUND
        elseif error isa FailedToProvePRLPFeasibility
            @debug "Failed to prove PRLP Feasibility"
            sepa.termination_status = FAILED_TO_PROVE_PRLP_FEASIBILITY
        elseif error isa FailedDisjunctiveLowerBoundTest
            @debug "Failed Disjunctive Lower Bound Test"
            sepa.termination_status = FAILED_DISJUNCTIVE_LOWER_BOUND_TEST
        elseif error isa PStarNotTight
            @debug "PStar Not Tight"
            sepa.termination_status = FAILED_TO_TIGHTEN_PSTAR
        elseif error isa BasestatZeroEncountered
            @debug "Basestat Zero Encountered"
            sepa.termination_status = BASESTAT_ZERO_ENCOUNTERED
        else
            rethrow(error)
        end

        # In any case if an error occurs, we return didnotrun to avoid future calls 
        return SCIP.SCIP_DIDNOTFIND
    end

    # Ordinary termination, return separated if cuts are found and didnotfind otherwise
    if sepa.termination_status == FOUND_CUTS
        return SCIP.SCIP_SEPARATED
    else
        return SCIP.SCIP_DIDNOTFIND
    end
end

function vpolyhedralcut_separation(sepa::VPCSeparator)
    # Aliasing for easier call
    scip = sepa.scipd
    start_time = time()

    # If cut limit is -1 or -2 convert them to the actual limit 
    # We do the conversion here because for option -2 we need the number of fractional variables
    if sepa.parameters.cut_limit == -1
        sepa.parameters.cut_limit = typemax(Int)
    elseif sepa.parameters.cut_limit == -2
        sepa.parameters.cut_limit = SCIP.SCIPgetNLPBranchCands(scip)
    end

    # Capture fractional variables statistic
    sepa.statistics.num_fractional_variables = SCIP.SCIPgetNLPBranchCands(scip)

    # Step 1: Construct NonBasicSpace and get LP Objective 
    lp_obj = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    nonbasic_space = NonBasicSpace(scip)
    # Capture LP Objective statistic
    sepa.statistics.lp_objective = lp_obj
    sepa.statistics.num_lp_rows = Int64(SCIP.SCIPgetNLPRows(scip))
    sepa.statistics.num_lp_columns = Int64(SCIP.SCIPgetNLPCols(scip))

    # Step 2: Get Disjunction
    @info "Getting Disjunction by Branch and Bound"
    lp_iter = Ref{Int64}()
    disjunction_timed = @timed get_disjunction_by_branchandbound(
        scip, sepa.parameters.n_leaves;
        log_path = joinpath(sepa.parameters.log_directory, "branch_and_bound.log"),
        time_limit = 0.5 * sepa.parameters.time_limit,
        lp_iter
    )
    disjunction = disjunction_timed.value
    sepa.statistics.branch_and_bound_lp_iteration = lp_iter[]
    sepa.statistics.branch_and_bound_time = disjunction_timed.time
    @info "Disjunction Obtained"

    # Step 3: Collect Point Ray
    @info "Collecting Points and Rays"
    point_ray_collection_timed = @timed get_point_ray_collection(
        scip, disjunction, nonbasic_space;
        time_limit = sepa.parameters.time_limit, start_time = start_time,
        lp_iter = lp_iter
    )
    point_ray_collection = point_ray_collection_timed.value
    sepa.statistics.point_ray_collection_time = point_ray_collection_timed.time
    sepa.statistics.point_ray_collection_lp_iterations = lp_iter[]
    sepa.statistics.num_points = length(get_points(point_ray_collection))
    sepa.statistics.num_rays = length(get_rays(point_ray_collection))
    @info "Finished Collecting Points and Rays"

    # Step 4: Test disjunctive_lower_bound
    # This test is always runned! The parameter sepa.parameters.test_disjunctive_lower_bound 
    # only determine whether we can skipped if the test fail 
    pstar = argmin(x -> get_objective_value(x), get_points(point_ray_collection))

    disjunctive_lower_bound = get_objective_value(pstar)
    sepa.statistics.disjunctive_lower_bound = disjunctive_lower_bound
    if is_LE(disjunctive_lower_bound, lp_obj) &&
        sepa.parameters.test_disjunctive_lower_bound
        @debug "Disjunctive Lower Bound is not strictly larger than the current LP"
        throw(FailedDisjunctiveLowerBoundTest())
    end

    # Step 5: Construct PRLP problem
    @info "Constructing PRLP"
    prlp_timed = @timed construct_prlp(
        scip,
        point_ray_collection;
        prlp_solve_method = PRLPsolveAlgorithm(sepa.parameters.prlp_solve_method),
        prlp_allow_warm_start = sepa.parameters.prlp_allow_warm_start
    )
    prlp = prlp_timed.value
    sepa.statistics.prlp_construction_time = prlp_timed.time
    sepa.statistics.prlp_num_columns = prlp.dimension
    sepa.statistics.prlp_num_rows = length(prlp.rays) + length(prlp.points)
    @info "PRLP Constructed"

    # Step 6: Gather separating solutions
    @info "Gathering Separating Solutions"
    status = Ref{GSSeturnCode}(GSS_OKAY)
    separating_solutions_timed = @timed gather_separating_solutions(
        prlp, point_ray_collection, status;
        cut_limit = sepa.parameters.cut_limit,
        time_limit = sepa.parameters.time_limit,
        start_time = start_time
    )
    if status[] == GSS_TIME_LIMIT
        sepa.termination_status = TIME_LIMIT_EXCEEDED_PRLP
    elseif status[] == GSS_CONSECUTIVE_FAIL
        sepa.termination_status = CONSECUTIVE_FAIL_LIMIT_REACHED
    end

    separating_solutions = separating_solutions_timed.value
    sepa.statistics.num_cuts = length(separating_solutions)
    sepa.statistics.prlp_separation_time = separating_solutions_timed.time
    sepa.statistics.objective_tried = length(prlp.solve_statistics)
    sepa.statistics.num_basis_restart = prlp.n_basis_restart
    @info "Basis Restart Count $(prlp.n_basis_restart)"
    # Capture statistics from PRLP
    sepa.statistics.prlp_solves_data = prlp.solve_statistics

    # Step 7: Process Cuts and Add to SCIP
    cutpool = CutPool(; variable_pointers = nonbasic_space.variable_pointers)
    for solution in separating_solutions
        cut = get_cut_from_separating_solution(
            solution, nonbasic_space
        )
        push!(cutpool, cut)
    end
    add_all_cuts!(scip, cutpool, sepa)

    # Sepa termination status have been modified due to Timelimit/ consecutive_fail
    if sepa.termination_status != NOT_RUN
        return nothing
    end

    if length(cutpool) > 0
        sepa.termination_status = FOUND_CUTS
    else
        sepa.termination_status = NO_CUTS_FOUND
    end
end