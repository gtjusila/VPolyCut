#
# structures.jl
# Structures needed for experiment
#

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