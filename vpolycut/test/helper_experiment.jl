using SCIP

"""
An experiment configuration object to store settings for an experiment. 

# Field 
- seperator::Bool Turn on separator? Default: True 
- heuristics::Bool Turn on heuristics? Default: True
- presolve::Bool Turn on presolve? Default: True
- propagation::Bool Turn on propagation? Default: True
- conflict::Bool Turn on conflict handling? Default: True
- zero_cut::Bool Allow cut with zero power? Default: True 
- symmetry::Bool Turn on symmetry handling? Default: True
- debug::Bool use SCIP debug solution? Default: False
- verbosity::Int Display verbosity level. Default: 5
- vpolycut::Bool should V Polyhedral Cut be used? Default: True
- node_limit::Int SCIP node limit. Default: -1
"""
@kwdef struct ExperimentConfiguration
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
    node_limit::Int64 = -1
end

"""
Solve `instance` with the given ExperimentConfiguration.
Instance should be located on the folder test\\data\\instances
Instance should be name instance.mps and if debug is turned on solution
should be named instance.sol
"""
function run_experiment(instance::String, config::ExperimentConfiguration)
    # STEP 1: Create a new model, Setup parameter and read problem
    optimizer = SCIP.Optimizer()
    inner = optimizer.inner

    # STEP 2: Set parameter
    
    setter = (par, val) -> SCIP.set_parameter(inner, par, val)
    turn_off_scip_heuristics(setter)
    
    if !config.separator turn_off_scip_separators(setter) end
    if !config.heuristics turn_off_scip_heuristics(setter) end
    if !config.presolve setter("presolving/maxrounds", 0) end
    if !config.propagation setter("propagating/maxroundsroot",0) end
    if !config.conflict setter("conflict/enable",false) end
    if config.zero_cut allow_zero_power_cut(setter) end
    if !config.symmetry setter("misc/usesymmetry", 0) end
    setter("display/verblevel",config.verbosity) 
    setter("limits/nodes",config.node_limit)

    # STEP 3: Load problem with debug solution (if debug is turned on)
    path = joinpath(@__DIR__,"data","instances",instance)
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(inner, path*".mps", C_NULL)
    if config.debug
        SCIP.set_parameter(inner,"misc/debugsol",joinpath(@__DIR__,"data","instances",instance*".sol"))
        SCIP.SCIPenableDebugSol(inner)
    end

    # STEP 4: Include V-Polyhedral Cut
    if config.vpolycut
        sepa = IntersectionSeparator(scipd = inner)
        SCIP.include_sepa(inner.scip[], inner.sepas, sepa; freq= 0)
    end

    SCIP.@SCIP_CALL SCIP.SCIPsolve(inner)
end

"""
Settings from SCIP
"""
function turn_off_scip_separators(setter::Function)
    setter("separating/disjunctive/freq", -1)
    setter("separating/impliedbounds/freq", -1)
    setter("separating/gomory/freq", -1)
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

function turn_off_scip_heuristics(setter::Function)
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

function allow_zero_power_cut(setter::Function)
    setter("separating/minefficacy", 0)
    setter("separating/minefficacyroot", 0)
    setter("cutselection/hybrid/minortho", 0)
    setter("cutselection/hybrid/minorthoroot", 0)
    setter("separating/poolfreq", 1)
end