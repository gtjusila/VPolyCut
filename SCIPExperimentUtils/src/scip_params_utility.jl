#
# setscipparam.jl
# Code to set scip parameter and turn on separators
#
using SCIP
using JuMP

export set_heuristics_emphasis_off, set_separators_emphasis_off, set_cut_selection_off, set_presolving_emphasis_off, set_root_node_propagation_off, set_everything_off, set_mode_cutting_plane_experiment

function get_parameter_setter_function(model::JuMP.AbstractModel)
    return (par, val) -> JuMP.set_attribute(model, par, val)
end

function set_heuristics_emphasis_off(model::JuMP.AbstractModel)
    setter = get_parameter_setter_function(model)

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
end

function set_separators_emphasis_off(model::JuMP.AbstractModel)
    setter = get_parameter_setter_function(model)

    # Turn off seperators
    setter("separating/disjunctive/freq", -1)
    setter("separating/impliedbounds/freq", -1)
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
end

function set_cut_selection_off(model::JuMP.AbstractModel)
    setter = get_parameter_setter_function(model)

    # Allow cut with zero efficacy
    setter("separating/minefficacy", 0)
    setter("separating/minefficacyroot", 0)
    setter("cutselection/hybrid/minortho", 0)
    setter("cutselection/hybrid/minorthoroot", 0)
    setter("cutselection/hybrid/minorthoroot", 0)
    setter("separating/filtercutpoolrel", false)
end

function set_presolving_emphasis_off(model::JuMP.AbstractModel)
    setter = get_parameter_setter_function(model)
    setter("presolving/maxrounds", 0)
end

function set_root_node_propagation_off(model::JuMP.AbstractModel)
    setter = get_parameter_setter_function(model)
    setter("propagating/maxroundsroot", 0)
end

function set_strong_branching_lookahead_off(model::JuMP.AbstractModel)
    setter = get_parameter_setter_function(model)
    setter("branching/relpscost/initcand", 0)
end
function set_everything_off(model::JuMP.AbstractModel)
    set_heuristics_emphasis_off(model)
    set_separators_emphasis_off(model)
    set_cut_selection_off(model)
    set_presolving_emphasis_off(model)
    set_root_node_propagation_off(model)
    set_strong_branching_lookahead_off(model)
end

function set_mode_cutting_plane_experiment(model::JuMP.AbstractModel)
    set_everything_off(model)

    setter = get_parameter_setter_function(model)

    # Set time limit
    setter("limits/time", 3600)

    # Only 1 round of separation
    setter("separating/maxroundsroot", 1)

    # Other settings
    setter("display/verblevel", 5)
    setter("limits/nodes", 1)
end