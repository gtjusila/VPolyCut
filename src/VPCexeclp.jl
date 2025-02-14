#
# This file contains the exec_lp function for VPC separation
#
using SCIP
using JuMP
using LinearAlgebra
import MathOptInterface as MOI
using Profile
using StatProfilerHTML
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
    @debug "VPC Separator Called"
    sepa.statistics.called += 1
    sepa.start_time = time()

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
    @debug "Finished checking necessary Preconditions"

    # Do everything in a try block to ensure time limit requirement
    @debug "Starting separation subroutine"
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
        if error isa TimeLimitExceeded
            @debug "Time Limit Exceeded"
            sepa.termination_status = TIME_LIMIT_EXCEEDED
        elseif error isa FailedToProvePRLPFeasibility
            @debug "Failed to prove PRLP Feasibility"
            sepa.termination_status = FAILED_TO_PROVE_PRLP_FEASIBILITY
        elseif error isa FailedDisjunctiveLowerBoundTest
            @debug "Failed Disjunctive Lower Bound Test"
            sepa.termination_status = FAILED_DISJUNCTIVE_LOWER_BOUND_TEST
        elseif error isa PStarNotTight
            @debug "PStar Not Tight"
            sepa.termination_status = FAILED_TO_TIGHTEN_PSTAR
        elseif error isa AssumptionViolated
            @debug "Assumption Violated"
            sepa.termination_status = ASSUMPTION_VIOLATED
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

    # If cut limit is -1 or -2 convert them to the actual limit 
    # We do the conversion here because for option -2 we need the number of fractional variables
    if sepa.parameters.cut_limit == -1
        sepa.parameters.cut_limit = typemax(Int)
    elseif sepa.parameters.cut_limit == -2
        sepa.parameters.cut_limit = SCIP.SCIPgetNLPBranchCands(scip)
    end
    # Capture fractional variables statistic
    sepa.statistics.n_fractional_variables = SCIP.SCIPgetNLPBranchCands(scip)

    # Step 1: Get LP objective Tableau
    sepa.lp_obj = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    sepa.tableau = construct_tableau_with_constraint_matrix(scip)
    sepa.nonbasic_space = NonBasicSpace(scip)
    sepa.lp_sol = get_solution_vector(scip)
    sepa.constraint_matrix = ConstraintMatrix(scip)

    # Step 2: Get Disjunction
    @debug "Getting Disjunction by Branch and Bound"
    sepa.disjunction = get_disjunction_by_branchandbound(scip, sepa.parameters.n_leaves)

    # Step 3: Collect Point Ray
    point_ray_time = time()
    sepa.point_ray_collection = get_point_ray_collection(
        scip, sepa.disjunction, sepa.nonbasic_space
    )
    @info "Point Ray Collection Time: $(time() - point_ray_time)"
    @debug "Number of points: $(num_points(sepa.point_ray_collection))"
    @debug "Number of rays: $(num_rays(sepa.point_ray_collection))"

    # Step 4: Test disjunctive_lower_bound
    # This test is always runned! The parameter sepa.parameters.test_disjunctive_lower_bound 
    # only determine whether we can skipped if the test fail 
    pstar = argmin(x -> get_objective_value(x), get_points(sepa.point_ray_collection))
    sepa.disjunctive_lower_bound = get_objective_value(pstar)
    if is_LE(sepa.disjunctive_lower_bound, sepa.lp_obj) &&
        sepa.parameters.test_disjunctive_lower_bound
        @debug "Disjunctive Lower Bound is not strictly larger than the current LP"
        throw(FailedDisjunctiveLowerBoundTest())
    end

    # Step 5: Construct PRLP problem
    prlp = construct_prlp(sepa.point_ray_collection)
    # Determine the method to solve PRLP
    @info "Checking Feasibility"
    if !PRLPcalibrate(prlp)
        throw(FailedToProvePRLPFeasibility())
    end

    # Step 6: Gather separating solutions
    separating_solution_time = time()
    separating_solutions = gather_separating_solutions(
        prlp, sepa.point_ray_collection;
        cut_limit = sepa.parameters.cut_limit,
        time_limit = sepa.parameters.time_limit,
        start_time = sepa.start_time
    )
    @info "Separating Solution Time: $(time() - separating_solution_time)"

    # Capture statistics from PRLP
    sepa.statistics.prlp_solves_data = prlp.solve_statistics

    # Step 7: Process Cuts and Add to SCIP
    sepa.cutpool = CutPool(; tableau = sepa.tableau)
    for solution in separating_solutions
        cut = get_cut_from_separating_solution(
            solution, sepa.constraint_matrix, sepa.nonbasic_space
        )
        push!(sepa.cutpool, cut)
    end
    add_all_cuts!(scip, sepa.cutpool, sepa)

    if length(sepa.cutpool) > 0
        sepa.termination_status = FOUND_CUTS
    else
        sepa.termination_status = NO_CUTS_FOUND
    end
end