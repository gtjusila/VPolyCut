import MathOptInterface as MOI
import Polyhedra

function hrep_to_constraint(hpoly,model,x)
    # Add Halfspaces
    for i in Polyhedra.halfspaces(hpoly)
        MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction(
                MOI.ScalarAffineTerm.(i.a,x),
                0.0,
            ),
            MOI.LessThan(i.β),
        )
    end
   # Add Hyperplanes
    for i in Polyhedra.hyperplanes(hpoly)
        MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction(
                MOI.ScalarAffineTerm.(i.a,x),
                0.0,
            ),
            MOI.EqualTo(i.β),
        )
    end 
end

function sepa_set_scip_parameters(setter::Function)
    
    # Disable SCIP Magic For Testing
    setter("display/verblevel", 5)

    setter("presolving/maxrounds", 0)
    setter("misc/usesymmetry", 0)

    setter("propagating/maxroundsroot",0)
    setter("propagating/maxrounds",0) 
    setter("conflict/enable",false)
    setter("constraints/linear/propfreq",-1)
    setter("constraints/linear/tightenboundsfreq",-1)
    setter("limits/nodes",1) #Only Do root nodes

    setter("separating/minefficacy", 0)
    setter("separating/minefficacyroot", 0)
    setter("cutselection/hybrid/minortho", 0)
    setter("cutselection/hybrid/minorthoroot", 0)
    setter("separating/poolfreq", 1)

    # These parameters come from `SCIP> set sepa emph off`
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
    setter("constraints/indicator/sepafreq", -1)

    # These parameters come from `SCIP> set heuristics off`
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