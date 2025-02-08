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

function SCIP.exec_lp(sepa::VPCSeparator)
    # Wrapper for writing logs
    logger = nothing
    if sepa.parameters.write_log
        logger = setup_file_logger(
            joinpath(sepa.parameters.log_directory, "vpc_separation.log")
        )
    else
        logger = ConsoleLogger(Logging.Debug)
    end

    with_logger(logger) do
        return _exec_lp(sepa)
    end
end

# Actual function
function _exec_lp(sepa::VPCSeparator)
    # Aliasing for easier call
    scip = sepa.scipd

    # Initialization: start timer, set separated to false, and increment call count 
    @debug "VPC Separator Called"
    sepa.statistics[CALLED] += 1
    sepa.start_time = time()
    sepa.separated = false

    # Check Preconditions and handle accordingly
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
    error_occurred = false
    try
        # Call Separation Subroutine
        vpolyhedralcut_separation(sepa)

        if sepa.separated
            sepa.termination_message = "CUT_FOUND"
            return SCIP.SCIP_SEPARATED
        else
            return SCIP.SCIP_DIDNOTFIND
        end
    catch error
        error_occurred = false

        if is_true(SCIP.SCIPinProbing(scip))
            SCIP.SCIPendProbing(scip)
        end

        if error isa TimeLimitExceeded
            @debug "Time Limit Exceeded"
            error_occurred = true
            sepa.termination_message = "TIME_LIMIT_EXCEEDED"
        elseif error isa FailedToProvePRLPFeasibility
            @debug "Failed to prove PRLP Feasibility"
            error_occurred = true
            sepa.termination_message = "FAILED_TO_PROVE_PRLP_FEASIBILITY"
        elseif error isa FailedDisjunctiveLowerBoundTest
            @debug "Failed Disjunctive Lower Bound Test"
            error_occurred = true
            sepa.termination_message = "FAILED_DISJUNCTIVE_LOWER_BOUND_TEST"
        elseif error isa PStarInfeasible
            @debug "PStar Is Infeasible"
            error_occurred = true
            sepa.termination_message = "PSTAR_INFEASIBLE"
        elseif error isa PStarNotTight
            @debug "PStar Not Tight"
            error_occurred = true
            sepa.termination_message = "PSTAR_NOT_TIGHT"
        else
            rethrow(error)
        end
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
    sepa.statistics[N_FRACTIONAL_VARIABLES] = SCIP.SCIPgetNLPBranchCands(scip)

    # Step 0: Get complemented tableau
    sepa.statistics[LP_OBJ] = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    construct_complemented_tableau(sepa)

    # Step 1: Get Disjunction
    @debug "Getting Disjunction by Branch and Bound"
    sepa.disjunction = get_disjunction_by_branchandbound(scip, sepa.parameters.n_leaves)

    #
    # Main Algorithm 
    #

    # Step 2: Setup projection used and get the points and rays from the disjunctions
    sepa.projection = create_projection_to_nonbasic_space(sepa.complemented_tableau)

    # Step 3: Collect Point Ray
    get_point_ray_collection(sepa)
    @debug "Number of points: $(num_points(sepa.point_ray_collection))"
    @debug "Number of rays: $(num_rays(sepa.point_ray_collection))"

    # Step 4: Setup cutpool
    @info "Tableau Density $(get_tableau_density(sepa.scipd, sepa.complemented_tableau.complemented_tableau))"
    sepa.cutpool = CutPool(; tableau=sepa.complemented_tableau, scip=scip)

    # Step 5: Verify that the disjunctive lower bound is strictly larger than the current LP
    sepa.disjunctive_lower_bound = minimum(
        x -> get_orig_objective_value(x), get_points(sepa.point_ray_collection)
    )
    if is_LE(scip, sepa.disjunctive_lower_bound, sepa.lp_obj)
        @debug "Disjunctive Lower Bound is not strictly larger than the current LP"
        throw(FailedDisjunctiveLowerBoundTest())
    end

    # Step 6: Solve Separation Problem
    solve_separation_subproblems(sepa)
end

function construct_complemented_tableau(sepa::VPCSeparator)
    original_tableau = construct_tableau_with_constraint_matrix(sepa.scipd)
    sepa.complemented_tableau = ComplementedTableau(original_tableau)
end
