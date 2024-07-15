import JuMP
import SCIP
import MathOptInterface as MOI
import VPolyCut

function setup_jump_model()
    model = setup_nocache_jump_model()
    return model
end

function include_separator(scip::SCIP.SCIPData, seperator::Type{T}; kwargs...) where {T<:SCIP.AbstractSeparator}
    sepa = seperator(; kwargs...)
    SCIP.include_sepa(scip.scip[], scip.sepas, sepa; freq=0, usessubscip=true)
end

"""
setup_jump_model

A JuMP model without cache. In standard JuMP model, we cannot control
when is the SCIP Optimizer object recreated so it will be impossible to do 
SCIP calls on it. 
"""
function setup_nocache_jump_model()
    inner = MOI.Bridges.full_bridge_optimizer(SCIP.Optimizer(), Float64)
    model = JuMP.direct_generic_model(Float64, inner)
    return model
end

function get_scipdata_from_model(model::T) where {T<:JuMP.AbstractModel}
    backend = JuMP.unsafe_backend(model)
    return backend.inner
end

function set_scip_parameters(model::T) where {T<:JuMP.AbstractModel}
    setter = (par, var) -> JuMP.set_attribute(model, par, var)

    # Turn off heuristics
    setter("heuristics/padm/freq", -1)
    setter("heuristics/ofins/freq", -1)
    setter("heuristics/trivialnegation/freq", -1)
    setter("heuristics/reoptsols/freq", -1)
    setter("heuristics/trivial/freq", -1)
    setter("heuristics/clique/freq", -1)
    setter("heuristics/locks/freq", -1)
    setter("heuristics/vbounds/freq", -1)
    setter("heuristics/shiftandpropagate/freq", -1)
    setter("heuristics/completesol/freq", -1)
    setter("heuristics/simplerounding/freq", -1)
    setter("heuristics/randrounding/freq", -1)
    setter("heuristics/zirounding/freq", -1)
    setter("heuristics/rounding/freq", -1)
    setter("heuristics/shifting/freq", -1)
    setter("heuristics/intshifting/freq", -1)
    setter("heuristics/oneopt/freq", -1)
    setter("heuristics/indicator/freq", -1)
    setter("heuristics/adaptivediving/freq", -1)
    setter("heuristics/farkasdiving/freq", -1)
    setter("heuristics/feaspump/freq", -1)
    setter("heuristics/conflictdiving/freq", -1)
    setter("heuristics/pscostdiving/freq", -1)
    setter("heuristics/fracdiving/freq", -1)
    setter("heuristics/nlpdiving/freq", -1)
    setter("heuristics/veclendiving/freq", -1)
    setter("heuristics/distributiondiving/freq", -1)
    setter("heuristics/objpscostdiving/freq", -1)
    setter("heuristics/rootsoldiving/freq", -1)
    setter("heuristics/linesearchdiving/freq", -1)
    setter("heuristics/guideddiving/freq", -1)
    setter("heuristics/rens/freq", -1)
    setter("heuristics/alns/freq", -1)
    setter("heuristics/rins/freq", -1)
    setter("heuristics/gins/freq", -1)
    setter("heuristics/lpface/freq", -1)
    setter("heuristics/crossover/freq", -1)
    setter("heuristics/undercover/freq", -1)
    setter("heuristics/subnlp/freq", -1)
    setter("heuristics/mpec/freq", -1)
    setter("heuristics/multistart/freq", -1)
    setter("heuristics/trysol/freq", -1)

    # Turn off seperators
    setter("separating/disjunctive/freq", -1)
    setter("separating/impliedbounds/freq", -1)
    setter("separating/gomory/freq", -1)
    #setter("separating/gmi/freq", -1)
    setter("separating/strongcg/freq", -1)
    setter("separating/aggregation/freq", -1)
    setter("separating/clique/freq", -1)
    setter("separating/zerohalf/freq", -1)
    setter("separating/mcf/freq", -1)
    setter("separating/flowcover/freq", -1)
    setter("separating/cmir/freq", -1)
    setter("separating/rapidlearning/freq", -1)
    setter("constraints/cardinality/sepafreq", -1)
    setter("constraints/SOS1/sepafreq", -1)
    setter("constraints/SOS2/sepafreq", -1)
    setter("constraints/varbound/sepafreq", -1)
    setter("constraints/knapsack/sepafreq", -1)
    setter("constraints/setppc/sepafreq", -1)
    setter("constraints/linking/sepafreq", -1)
    setter("constraints/or/sepafreq", -1)
    setter("constraints/and/sepafreq", -1)
    setter("constraints/xor/sepafreq", -1)
    setter("constraints/linear/sepafreq", -1)
    setter("constraints/orbisack/sepafreq", -1)
    setter("constraints/symresack/sepafreq", -1)
    setter("constraints/logicor/sepafreq", -1)
    setter("constraints/cumulative/sepafreq", -1)
    setter("constraints/nonlinear/sepafreq", -1)
    setter("separating/mixing/freq", -1)
    setter("separating/rlt/freq", -1)
    setter("constraints/indicator/sepafreq", -1)

    # Allow cut with zero efficacy
    setter("separating/minefficacy", 0)
    setter("separating/minefficacyroot", 0)
    setter("cutselection/hybrid/minortho", 0)
    setter("cutselection/hybrid/minorthoroot", 0)
    setter("cutselection/hybrid/minorthoroot", 0)
    setter("separating/filtercutpoolrel", false)

    # For stability disable restart
    setter("limits/restarts", 0)

    # Set time limit
    setter("limits/time", 3600)

    # Only 1 round of separation
    setter("separating/maxroundsroot", 1)

    # Other settings
    setter("display/verblevel", 5)
    setter("limits/nodes", 1)
    setter("branching/relpscost/initcand", 0)
end

function set_scip_parameters_easy(model::T) where {T<:JuMP.AbstractModel}
    set_scip_parameters(model::T)
    JuMP.set_attribute(model, "propagating/maxroundsroot", 0)
    JuMP.set_attribute(model, "presolving/maxrounds", 0)
end