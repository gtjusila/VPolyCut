using Test
import Polyhedra
import SCIP
import MathOptInterface as MOI

include("./utils.jl")
include("./intersection.jl")

function run_test_miplib() 

    # STEP 1: Create a new model, Setup parameter and read problem
    optimizer = SCIP.Optimizer()
    inner = optimizer.inner
    setter = (par, val) -> SCIP.set_parameter(inner, par, val)
    turn_off_scip_miscellaneous(setter)
    turn_off_scip_heuristics(setter)
    turn_off_scip_separators(setter)
    allow_zero_power_cut(setter)
    setter("display/verblevel",5)
    setter("limits/nodes",1)
    setter("separating/maxstallroundsroot",1000) 

    SCIP.@SCIP_CALL SCIP.SCIPreadProb(inner.scip[],joinpath(@__DIR__,"data","gen-ip054.mps"),C_NULL)
    SCIP.set_parameter(inner,"misc/debugsol",joinpath(@__DIR__,"data","gen-ip054.sol"))
    SCIP.SCIPenableDebugSol(inner.scip[]) 

    # Add our separator to scip
    sepa = IntersectionSeparator(scipd = inner, debug_sol_path = joinpath(@__DIR__,"data","gen-ip054.sol"))
    SCIP.include_sepa(inner.scip[], inner.sepas, sepa; freq= 0) 
    
    # Do the actual solve
    SCIP.@SCIP_CALL SCIP.SCIPsolve(inner.scip[])

end

run_test_miplib()