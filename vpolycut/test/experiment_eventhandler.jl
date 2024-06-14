using SCIP

include("experiment_structures.jl")

@kwdef mutable struct RecordSolutionAfterFirstRootLPSolve <: SCIP.AbstractEventHandler
    scipd::SCIP.SCIPData
    experiment_log::IO
    store::ExperimentStore
end

function SCIP.eventexec(event::RecordSolutionAfterFirstRootLPSolve)
    scip = event.scipd
    store = event.store

    if (SCIP.SCIPgetCurrentNode(scip) == SCIP.SCIPgetRootNode(scip)) 
        @info "First LP Solved" 
    end

    sol = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateLPSol(scip, sol, C_NULL)
    
    objval = SCIP.SCIPgetSolOrigObj(scip, sol[])

    store.firstlpobj = objval
    store.initialgap = abs(store.refobj-store.firstlpobj)
    
    nfrac = SCIP.SCIPgetNLPBranchCands(scip)
    frac = Ref{Ptr{SCIP.SCIP_Real}}(C_NULL) 
    
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBranchCands(scip, C_NULL, C_NULL, frac,C_NULL, C_NULL, C_NULL)
    frac = unsafe_wrap(Vector{SCIP.SCIP_Real},frac[],nfrac)
    
    @warn nfrac
    println(frac)
    
    store.nfrac = nfrac
    SCIP.@SCIP_CALL SCIP.SCIPfreeSol(scip, sol)
end

"""
Add RecordSolutionAfterFirstRootLPSolve event handler to scip
"""
function add_first_lp_event_handler(scip::SCIP.SCIPData, experiment_log::IO, store::ExperimentStore)
    firstlp = RecordSolutionAfterFirstRootLPSolve(scipd = scip, experiment_log = experiment_log, store = store)
    SCIP.include_event_handler(scip.scip[], scip.eventhdlrs, firstlp)
    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(scip)
    SCIP.@SCIP_CALL SCIP.SCIPcatchEvent(
        scip.scip[], 
        SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, 
        scip.eventhdlrs[firstlp],
        C_NULL,
        C_NULL
    ) 
end