using SCIP
using VPolyCut

@kwdef mutable struct TableauTest <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
    callback::Function
end

function SCIP.exec_lp(sepa::TableauTest)
    sepa.callback(sepa.scipd)
    return SCIP.SCIP_DIDNOTRUN
end
