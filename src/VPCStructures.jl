using SCIP
using JuMP
using MathOptInterface

"""
VPolyhedral Cut Separator

Implementation of Algorithm 1 from 
Balas, Egon, and Aleksandr M. Kazachkov. "V-polyhedral disjunctive cuts." arXiv preprint arXiv:2207.13619 (2022).
Disjunction are obtained from partial branch and bound trees
"""

"""
VPCParameters

A struct to store the parameters of the VPCSeparator. An instance of this class
is passed during the creation of the VPCSeparator.
"""
@kwdef mutable struct VPCParameters
    "Maximum number of separation rounds"
    call_limit::Int = 1
    "Maximum number of cuts to generate in a round. -1 for no limit, -2 for the number of fractional variable in LP Solution"
    cut_limit::Int = -2
    "Directory to write cut log"
    log_directory::String = ""
    "Max rounds to participate in"
    max_round::Int = 1
    "Number of leaves in the disjunction"
    n_leaves::Int = 2
    "Apply Beta Scaling"
    apply_beta_scaling::Bool = true
    "PRLP allow warm start"
    prlp_allow_warm_start::Bool = true
    "PRLP solve method"
    prlp_solve_method::Int = 1
    "PRLP maximum consecutive fail"
    prlp_max_consecutive_fail::Int = 5
    "PRLP minimum increase to be non stagnating"
    prlp_min_increase_non_stagnating::Float64 = 0.0
    "PRLP maximum number of stagnating objective"
    prlp_max_consecutive_stagnation::Int = 10

    "Test if disjunctive_lower_bound is better than LP Objective"
    test_disjunctive_lower_bound::Bool = true
    "Time Limit"
    time_limit::Float64 = typemax(64)
end

@kwdef mutable struct VPCStatistics
    branch_and_bound_lp_iterations::Int64 = 0.0
    branch_and_bound_time::Float64 = 0.0

    called::Int = 0
    cbar_test::Bool = true
    disjunctive_lower_bound::SCIP.SCIP_Real = 0.0
    disjunctive_lower_bound_history::Vector{SCIP.SCIP_Real} = []

    lp_objective::SCIP.SCIP_Real = 0.0

    num_basis_restart::Int = 0
    num_cuts::Int = 0
    num_fractional_variables::Int = 0
    num_lp_columns::Int = 0
    num_lp_rows::Int = 0
    num_points::Int = 0
    num_rays::Int = 0

    point_ray_collection_time::Float64 = 0.0
    point_ray_collection_lp_iterations::Int64 = 0

    prlp_construction_time::Float64 = 0.0
    prlp_num_columns::Int = 0
    prlp_num_rows::Int = 0
    prlp_separation_time::Float64 = 0.0
    prlp_solves_data::Vector{Any} = []
    prlp_solve_method::String = "PRIMAL_SIMPLEX"
    prlp_percent_disjunctive_gap_closed_history::Vector{SCIP.SCIP_Real} = []

    root_lp_iterations::Int64 = 0
    total_time_taken::Float64 = 0.0
    objective_tried::Int = 0
end

"""
    VPCTerminationStatus

An enum for possible termination status of the VPCSeparator
"""
@enum VPCTerminationStatus begin
    BASESTAT_ZERO_ENCOUNTERED
    CONSECUTIVE_FAIL_LIMIT_REACHED
    FAILED_DISJUNCTIVE_LOWER_BOUND_TEST
    FAILED_TO_PROVE_PRLP_FEASIBILITY
    FAILED_TO_TIGHTEN_PSTAR
    FOUND_CUTS
    LP_ERROR
    NO_CUTS_FOUND
    NOT_RUN
    TIME_LIMIT_EXCEEDED_BRANCHANDBOUND
    TIME_LIMIT_EXCEEDED_COLLECTION
    TIME_LIMIT_EXCEEDED_PRLP
end

@kwdef mutable struct VPCSharedData
    cutpool::Union{Nothing,CutPool} = nothing
    disjunction::Union{Nothing,Disjunction} = nothing
    disjunctive_lower_bound::SCIP.SCIP_Real = 0.0
    lp_obj::Float64 = 0.0
    lp_obj_nonbasic::Float64 = 0.0
    original_points::Vector{Point} = []
    nonbasic_space::Union{Nothing,NonBasicSpace} = nothing
    point_ray_collection::Union{Nothing,PointRayCollection} = nothing
    prlp::Union{Nothing,PRLP} = nothing
    scipd::SCIP.SCIPData
    separating_solutions::Union{Nothing,Vector{Vector{SCIP.SCIP_Real}}} = nothing
    start_time::Float64 = 0.0
end

"""
VPCSeparator

A class to hold the separator data. The actual separator object as according to the SCIP standard.

Constructors:
- `VPCSeparator(scipd::SCIP.SCIPData, params::VPCParameters)`: Create a new VPCSeparator
"""
@kwdef mutable struct VPCSeparator <: SCIP.AbstractSeparator
    # Shared Data among the functions
    "Shared VPC Data"
    shared_data::VPCSharedData
    "SEPA Parameters"
    parameters::VPCParameters
    "SEPA statistics"
    statistics::VPCStatistics = VPCStatistics()

    "should be skipped?"
    should_be_skipped::Bool = false

    # Return message
    "Termination Message"
    termination_status::VPCTerminationStatus = NOT_RUN
end

# Constructor
function VPCSeparator(params::VPCParameters)
    obj = VPCSeparator(; parameters = params)
    return obj
end

# Include Helper
function include_vpolyhedral_sepa(
    scipd::SCIP.SCIPData;
    parameters::VPCParameters = VPCParameters(),
    priority::Integer = 9999,
    freq::Integer = 1,
    maxbounddist::Real = 0.0,
    delay::Bool = false
)
    sepa = VPCSeparator(;
        shared_data = VPCSharedData(; scipd = scipd),
        parameters
    )
    SCIP.include_sepa(
        scipd.scip[], scipd.sepas, sepa; priority = priority,
        freq = freq, maxbounddist = maxbounddist,
        delay = delay, usessubscip = false
    )
    return sepa
end
