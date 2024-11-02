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
    "Should zeroing heuristic be used?"
    zeroing_heuristic::Bool = false
    "Should log be written?"
    write_log::Bool = false
    "Directory to write cut log"
    log_directory::String = ""
    "HiGHS LP Method"
    lp_solving_method::Int = 4
end

"""
VPCSeparator

A class to hold the separator data. The actual separator object as according to the SCIP standard.

Constructors:
- `VPCSeparator(scipd::SCIP.SCIPData, params::VPCParameters)`: Create a new VPCSeparator
"""
@kwdef mutable struct VPCSeparator <: SCIP.AbstractSeparator
    "Pointer to SCIP"
    scipd::SCIP.SCIPData
    "Number of times the seperation routine have been called"
    called::Int = 0
    "Have any cut been found during the round?"
    separated::Bool = false
    "SEPA Parameters"
    parameters::VPCParameters

    # Statistics
    "Termination Message"
    termination_message::String = ""
    "Disjunctive Lower Bound"
    disjunctive_lower_bound::SCIP.SCIP_Real = 0.0
    "Number of Fractional Variables"
    n_fractional_variables::Int = 0
    "PRLP Solves Statistics"
    prlp_solves::Vector{Dict{String,Any}} = []

    "Complemented Tableau"
    complemented_tableau::Union{Nothing,ComplementedTableau} = nothing
    "Disjunction"
    disjunction::Vector{Node} = []
    "PointRayCollection"
    point_ray_collection::Union{Nothing,PointRayCollection} = nothing
    "Projection"
    projection::Union{Nothing,Projection} = nothing
    "Cut Pool"
    cutpool::Union{Nothing,CutPool} = nothing
    "Separating Problem"
    separating_problem::Union{Nothing,AbstractModel} = nothing
    "LP Solution"
    lp_sol::Union{Nothing,Vector{SCIP.SCIP_Real}} = nothing
    "LP solution objective value"
    lp_obj::Union{Nothing,SCIP.SCIP_Real} = nothing
end

# Constructor
function VPCSeparator(scipd::SCIP.SCIPData, params::VPCParameters)
    return VPCSeparator(; scipd=scipd, parameters=params)
end

# Include Helper
function include_vpolyhedral_sepa(
    scipd::SCIP.SCIPData;
    n_leaves=2,
    cut_limit=-2,
    write_log=false,
    log_directory="",
    zeroing_heuristic=false,
    lp_solving_method=4
)
    parameters = VPCParameters(;
        n_leaves=n_leaves,
        cut_limit=cut_limit,
        write_log=write_log,
        log_directory=log_directory,
        zeroing_heuristic=zeroing_heuristic,
        lp_solving_method=lp_solving_method)

    if write_log
        if !isdir(log_directory)
            mkpath(log_directory)
        end
    end

    sepa = VPCSeparator(scipd, parameters)
    SCIP.include_sepa(
        scipd.scip[], scipd.sepas, sepa; priority=9999, freq=0, usessubscip=true
    )
    return sepa
end
