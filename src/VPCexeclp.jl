#
# This file contains the exec_lp function for VPC separation
#
using SCIP
using JuMP
using HiGHS
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

    # Increase the call counter and set separated to false (since no cut have been found)
    @debug "VPC Separator Called"
    sepa.called += 1
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
    try
        # Call Separation Subroutine
        vpolyhedralcut_separation(sepa)
    catch e
        print(e)
        if SCIP.SCIPinProbing(scip) == 1
            SCIP.SCIPendProbing(scip)
        end
        if e isa TimeLimitExceeded
            # Do nothing
        elseif e isa FailedToProvePRLPFeasibility
            @debug "Failed to prove PRLP Feasibility"
        else
            rethrow(e)
        end
    end

    if sepa.separated
        return SCIP.SCIP_SEPARATED
    else
        return SCIP.SCIP_DIDNOTRUN
    end
end

function vpolyhedralcut_separation(sepa::VPCSeparator)
    #
    # Preparation
    #
    scip = sepa.scipd

    # If cut limit is -1 or -2 convert them to the actual limit 
    # We do the conversion here because for option -2 we need the number of fractional variables
    if sepa.parameters.cut_limit == -1 || sepa.parameters.cut_limit == -2
        sepa.parameters.cut_limit = get_cut_limit(sepa, sepa.parameters)
    end

    # Step 0: Get complemented tableau
    sepa.lp_obj = SCIP.SCIPgetLPObjval(scip)
    construct_complemented_tableau(sepa)

    # Step 1: Get Disjunction
    @debug "Getting Disjunction by Branch and Bound"
    get_disjunction_by_branchandbound(sepa)

    #
    # Main Algorithm 
    #

    # Step 2: Setup projection used and get the points and rays from the disjunctions
    sepa.projection = create_projection_to_nonbasic_space(sepa.complemented_tableau)
    get_point_ray_collection(sepa)
    @debug "Number of points: $(num_points(sepa.point_ray_collection))"
    @debug "Number of rays: $(num_rays(sepa.point_ray_collection))"

    # Step 3: Setup Cut Pool 
    sepa.cutpool = CutPool(; tableau=sepa.complemented_tableau, scip=scip)

    # Step 4: Solve Separation Problem
    solve_separation_subproblems(sepa)
end

function get_cut_limit(
    sepa::VPCSeparator, sepa_parameter::VPCParameters
)
    if sepa_parameter.cut_limit == -1
        return typemax(Int)
    elseif sepa_parameter.cut_limit == -2
        return SCIP.SCIPgetNLPBranchCands(sepa.scipd)
    else
        return sepa_parameter.cut_limit
    end
end

function construct_complemented_tableau(sepa::VPCSeparator)
    original_tableau = construct_tableau_with_constraint_matrix(sepa.scipd)
    vars = map(1:get_nvars(original_tableau)) do i
        return get_var_from_column(original_tableau, i)
    end
    for var in vars
        if is_EQ(sepa.scipd, get_ub(var),get_lb(var))
            continue
        end
        if get_basis_status(var) == SCIP.SCIP_BASESTAT_UPPER
            @debug get_obj(var) get_basis_status(var) get_lb(var) get_ub(var) get_sol(var) get_obj(var)
            @assert is_LE(sepa.scipd,get_obj(var),0.0)
        end
        if get_basis_status(var) == SCIP.SCIP_BASESTAT_LOWER
            @debug get_obj(var) get_basis_status(var) get_lb(var) get_ub(var) get_sol(var) get_obj(var)
            @assert is_GE(sepa.scipd, get_obj(var), 0.0)
        end
    end
    @debug "LP Objective is $(dot(c_bar, lp_sol)) SepaLP OBJ is $(sepa.lp_obj)"
    sepa.complemented_tableau = ComplementedTableau(original_tableau)
    @debug "Complemented $(length(get_complemented_columns(sepa.complemented_tableau))) columns"
    vars = map(1:get_nvars(original_tableau)) do i
        return get_var_from_column(sepa.complemented_tableau, i)
    end
    for var in vars
        if is_EQ(sepa.scipd, get_ub(var),get_lb(var))
            continue
        end
        @assert get_basis_status(var) != SCIP.SCIP_BASESTAT_UPPER 
        if get_basis_status(var) == SCIP.SCIP_BASESTAT_LOWER
            @debug get_obj(var) get_basis_status(var) get_lb(var) get_ub(var) get_sol(var) get_obj(var)
            @assert is_GE(sepa.scipd, get_obj(var), 0.0)
        end
    end
    return sepa.complemented_tableau
end

function get_disjunction_by_branchandbound(sepa::VPCSeparator)
    branchandbound = BranchAndBound(
        sepa.scipd; max_leaves=sepa.parameters.n_leaves
    )
    execute_branchandbound(branchandbound)
    sepa.disjunction = get_leaves(branchandbound)
end