using SCIP

@kwdef mutable struct RecordSolutionAfterFirstRootLPSolve <: SCIP.AbstractEventHandler
    scipd::SCIP.SCIPData
end

function SCIP.eventexec(event::RecordSolutionAfterFirstRootLPSolve)
   if (SCIP.SCIPgetCurrentNode(event.scipd) == SCIP.SCIPgetRootNode(event.scipd)) 
        @warn "First LP Solved" 
   end
end