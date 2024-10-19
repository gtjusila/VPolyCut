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
    # Aliasing for easier call
    scip = sepa.scipd
    sepa.called += 1

    # Check Preconditions and handle accordingly
    if SCIP.SCIPgetStage(scip) != SCIP.SCIP_STAGE_SOLVING
        return SCIP.SCIP_DIDNOTRUN
    end
    if SCIP.SCIPisLPSolBasic(scip) == 0
        return SCIP.SCIP_DELAYED
    end
    if SCIP.SCIPgetLPSolstat(scip) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        return SCIP.SCIP_DELAYED
    end

    # Call Separation Routine
    vpolyhedralcut_separation(sepa)

    if sepa.separated
        return SCIP.SCIP_SEPARATED
    else
        return SCIP.SCIP_DIDNOTFIND
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
    construct_complemented_tableau(sepa)

    # Step 1: Get Disjunction
    get_disjunction_by_branchandbound(sepa)

    #
    # Main Algorithm 
    #

    # Step 2: Setup projection used and get the points and rays from the disjunctions
    sepa.projection = create_projection_to_nonbasic_space(sepa.complemented_tableau)
    get_point_ray_collection(sepa)
    with_logger(ConsoleLogger()) do
        @info "Number of points: $(num_points(sepa.point_ray_collection))"
        @info "Number of rays: $(num_rays(sepa.point_ray_collection))"
    end

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
    sepa.complemented_tableau = ComplementedTableau(original_tableau)
    return sepa.complemented_tableau
end

function get_disjunction_by_branchandbound(sepa::VPCSeparator)
    @info "Starting Branch And Bound to generate disjunction"
    branchandbound = BranchAndBound(
        sepa.scipd; max_leaves=sepa.parameters.n_leaves
    )
    execute_branchandbound(branchandbound)
    sepa.disjunction = get_leaves(branchandbound)
    @info "Branch and bound completed"
end