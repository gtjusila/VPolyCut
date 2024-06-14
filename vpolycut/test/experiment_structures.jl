"""
An experiment configuration object to store settings for an experiment. 

# Field 
- path::String The path to write any file output. The only required parameter.
- type::String The type of the experiment
- seperator::Bool Turn on separator? Default: True 
- heuristics::Bool Turn on heuristics? Default: True
- presolve::Bool Turn on presolve? Default: True
- propagation::Bool Turn on propagation? Default: True
- conflict::Bool Turn on conflict handling? Default: True
- zero_cut::Bool Allow cut with zero power? Default: True 
- symmetry::Bool Turn on symmetry handling? Default: True
- debug::Bool use debug solution? Default: False
- verbosity::Int Display verbosity level. Default: 5
- log_file::string Directory to write the log file. Default: Base.stdout 
- vpolycut::Bool should V Polyhedral Cut be used? Default: True
- vpolycut_limit::Int64 Maximum number of times the V Polyhedral Cut Seperator can be called? (-1 for no limit) Default: -1
- gomory::Bool Override seperator settings and turn on Gomory Cut for the root node? Default: False 
- node_limit::Int SCIP node limit. (-1 for no limit) Default: -1
- max_separounds::Int Maximum separator rounds in root Default: -1
"""
@kwdef struct ExperimentConfiguration
    path::String
    type::String = "default"
    separator::Bool = true
    heuristics::Bool = true
    presolve::Bool = true
    propagation::Bool = true 
    conflict::Bool = true 
    zero_cut::Bool = true
    symmetry::Bool = true
    debug::Bool = false
    verbosity::Int64 = 5
    vpolycut::Bool = true
    vpolycut_limit::Int64 = -1
    gomory::Bool = false
    node_limit::Int64 = -1
    max_separounds::Int64 = -1
end

"""
An ExperimentStore to record data from the experiment
"""
@kwdef mutable struct ExperimentStore
    refobj::SCIP.SCIP_Real = 0.0
    firstlpobj::SCIP.SCIP_Real = 0.0
    initialgap::SCIP.SCIP_Real = 0.0
    finalgap::SCIP.SCIP_Real = 0.0
    finallpobj::SCIP.SCIP_Real = 0.0
    gapclosed::SCIP.SCIP_Real = 0.0
    nfrac::Int64 = 0
end