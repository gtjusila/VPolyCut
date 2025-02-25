#
# This file contains the exec_lp function for VPC separation
#
using SCIP
using LinearAlgebra

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
    scip = sepa.shared_data.scipd

    # Initialization: start timer, set separated to false, and increment call count 
    @info "VPC Separator Called"
    sepa.statistics.called += 1

    # Check Preconditions and handle accordingly
    @info "Checking Precondition"
    if sepa.should_be_skipped
        return SCIP.SCIP_DIDNOTRUN
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

    # Setup and capture necessary statistic
    if !is_numeric_scip_set()
        set_numeric_scip(scip)
    end
    # If cut limit is -1 or -2 convert them to the actual limit 
    # We do the conversion here because for option -2 we need the number of fractional variables
    if sepa.parameters.cut_limit == -1
        sepa.parameters.cut_limit = typemax(Int)
    elseif sepa.parameters.cut_limit == -2
        sepa.parameters.cut_limit = SCIP.SCIPgetNLPBranchCands(scip)
    end
    # Capture fractional variables statistic and root node lp iterations
    sepa.statistics.num_fractional_variables = SCIP.SCIPgetNLPBranchCands(scip)
    sepa.statistics.root_lp_iterations = SCIP.SCIPgetNRootFirstLPIterations(scip)

    # Do everything in a try block 
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
        #sepa.should_be_skipped = true

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
            sepa.should_be_skipped = true
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
    scip = sepa.shared_data.scipd
    sepa.shared_data.start_time = time()

    # Step 1: Construct NonBasicSpace and get LP Objective 
    sepa.shared_data.lp_obj = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    sepa.shared_data.nonbasic_space = NonBasicSpace(scip)

    # Capture LP Objective statistic
    sepa.statistics.lp_objective = sepa.shared_data.lp_obj
    sepa.statistics.num_lp_rows = Int64(SCIP.SCIPgetNLPRows(scip))
    sepa.statistics.num_lp_columns = Int64(SCIP.SCIPgetNLPCols(scip))

    # Step 2: Get Disjunction
    @info "Getting Disjunction by Branch and Bound"
    disjunction_timed = @timed get_disjunction_by_branchandbound(sepa)
    sepa.shared_data.disjunction = disjunction_timed.value
    sepa.statistics.branch_and_bound_time = disjunction_timed.time
    @info "Disjunction Obtained"

    # Step 3: Collect Point Ray
    @info "Collecting Points and Rays"
    point_ray_collection_timed = @timed get_point_ray_collection(sepa)
    sepa.shared_data.point_ray_collection = point_ray_collection_timed.value
    sepa.statistics.point_ray_collection_time = point_ray_collection_timed.time
    sepa.statistics.num_points = length(get_points(sepa.shared_data.point_ray_collection))
    sepa.statistics.num_rays = length(get_rays(sepa.shared_data.point_ray_collection))
    @info "Finished Collecting Points and Rays"

    # Step 4: Test disjunctive_lower_bound
    # This test is always runned! The parameter sepa.parameters.test_disjunctive_lower_bound 
    # only determine whether we can skipped if the test fail 
    pstar = argmin(
        x -> get_objective_value(x), get_points(sepa.shared_data.point_ray_collection)
    )
    disjunctive_lower_bound = get_objective_value(pstar)
    sepa.statistics.disjunctive_lower_bound = disjunctive_lower_bound
    @info "Disjunctive Lower bound $(disjunctive_lower_bound)"
    if is_LE(disjunctive_lower_bound, sepa.shared_data.lp_obj) &&
        sepa.parameters.test_disjunctive_lower_bound
        @debug "Disjunctive Lower Bound is not strictly larger than the current LP"
        throw(FailedDisjunctiveLowerBoundTest())
    end

    # Step 5: Construct PRLP problem
    @info "Constructing PRLP"
    prlp_timed = @timed construct_prlp(sepa)
    sepa.shared_data.prlp = prlp_timed.value
    sepa.statistics.prlp_construction_time = prlp_timed.time
    sepa.statistics.prlp_num_columns = sepa.shared_data.prlp.dimension
    sepa.statistics.prlp_num_rows =
        length(sepa.shared_data.prlp.rays) + length(sepa.shared_data.prlp.points)
    @info "PRLP Constructed"

    # Step 6: Gather separating solutions
    @info "Gathering Separating Solutions"
    separating_solutions_timed = @timed gather_separating_solutions(sepa)

    sepa.shared_data.separating_solutions = separating_solutions_timed.value
    sepa.statistics.num_cuts = length(sepa.shared_data.separating_solutions)
    sepa.statistics.prlp_separation_time = separating_solutions_timed.time
    sepa.statistics.objective_tried = length(sepa.shared_data.prlp.solve_statistics)
    sepa.statistics.num_basis_restart = sepa.shared_data.prlp.n_basis_restart
    @info "Basis Restart Count $(sepa.shared_data.prlp.n_basis_restart)"
    # Capture statistics from PRLP
    sepa.statistics.prlp_solves_data = sepa.shared_data.prlp.solve_statistics

    # Step 7: Process Cuts and Add to SCIP
    sepa.shared_data.cutpool = CutPool(;
        variable_pointers = sepa.shared_data.nonbasic_space.variable_pointers
    )
    for solution in sepa.shared_data.separating_solutions
        cut = get_cut_from_separating_solution(
            solution, sepa.shared_data.nonbasic_space
        )
        push!(sepa.shared_data.cutpool, cut)
    end
    add_all_cuts!(scip, sepa.shared_data.cutpool, sepa)

    # Sepa termination status have been modified due to Timelimit/ consecutive_fail
    if sepa.termination_status != NOT_RUN
        return nothing
    end

    if length(sepa.shared_data.cutpool) > 0
        sepa.termination_status = FOUND_CUTS
    else
        sepa.termination_status = NO_CUTS_FOUND
    end
end