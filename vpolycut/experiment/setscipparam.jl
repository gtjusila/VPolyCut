#
# setscipparam.jl
# Code to set scip parameter and turn on separators
#
import SCIP

function setscipsettings(setter::Function)
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
    setter("separating/gmi/freq", -1)
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

    # Uncomment the following for super easy instances
    setter("propagating/maxroundsroot", 0)
    setter("presolving/maxrounds", 0)
end

function includegomorysepa(setter::Function)
    setter("separating/gmi/freq", 0)
    setter("separating/gmi/priority", 9999)
    setter("separating/gmi/maxsuppabs", 5000)
    setter("separating/gmi/dynamiccuts", false)
    setter("separating/gmi/maxsupprel", 1.0)
    setter("separating/gmi/forcecuts", true)
end

function includevpcsepa(scip::SCIP.SCIPData)
    sepa = IntersectionSeparator(scipd=scip)
    SCIP.include_sepa(scip.scip[], scip.sepas, sepa; freq=0, usessubscip=true)
end