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
    "Number of leaves in the disjunction"
    n_leaves::Int = 2
    "Maximum number of cuts to generate in a round. -1 for no limit, -2 for the number of fractional variable in LP Solution"
    cut_limit::Int = -2
    "Maximum number of separation rounds"
    call_limit::Int = 1
    "Directory to write cut log"
    log_directory::String = ""
    "Time Limit"
    time_limit::Float64 = typemax(64)
    "Test if disjunctive_lower_bound is better than LP Objective"
    test_disjunctive_lower_bound::Bool = true
end

@kwdef mutable struct VPCStatistics
    called::Int = 0
    prlp_solves_data::Vector{Any} = []
    n_fractional_variables::Int = 0
    cbar_test::Bool = true
    point_ray_collection_time::Float64 = 0.0
    n_points::Int = 0
    n_rays::Int = 0
    prlp_solve_method::String = "PRIMAL_SIMPLEX"
    disjunctive_lower_bound::Float64 = 0.0
    lp_objective::Float64 = 0.0
    prlp_separation_time::Float64 = 0.0
    prlp_construction_time::Float64 = 0.0
    prlp_feasibility_proving_time::Float64 = 0.0
    number_of_cuts::Int = 0
end

@enum VPCTerminationStatus begin
    TIME_LIMIT_EXCEEDED
    FAILED_TO_PROVE_PRLP_FEASIBILITY
    FAILED_DISJUNCTIVE_LOWER_BOUND_TEST
    BASESTAT_ZERO_ENCOUNTERED
    FAILED_TO_TIGHTEN_PSTAR
    ASSUMPTION_VIOLATED
    FOUND_CUTS
    NO_CUTS_FOUND
    NOT_RUN
end

"""
VPCSeparator

A class to hold the separator data. The actual separator object as according to the SCIP standard.

Constructors:
- `VPCSeparator(scipd::SCIP.SCIPData, params::VPCParameters)`: Create a new VPCSeparator
"""
@kwdef mutable struct VPCSeparator <: SCIP.AbstractSeparator
    # Shared Data among the functions
    "SCIP Data"
    scipd::SCIP.SCIPData
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
function VPCSeparator(scipd::SCIP.SCIPData, params::VPCParameters)
    obj = VPCSeparator(; scipd = scipd, parameters = params)
    return obj
end

# Include Helper
function include_vpolyhedral_sepa(
    scipd::SCIP.SCIPData;
    n_leaves = 2,
    cut_limit = -2,
    log_directory = "",
    time_limit = typemax(Float64)
)
    parameters = VPCParameters(;
        n_leaves = n_leaves,
        cut_limit = cut_limit,
        log_directory = log_directory,
        time_limit = time_limit
    )

    sepa = VPCSeparator(scipd, parameters)
    SCIP.include_sepa(
        scipd.scip[], scipd.sepas, sepa; priority = 9999, freq = 0, usessubscip = true
    )
    return sepa
end
