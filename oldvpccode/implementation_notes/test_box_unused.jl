using Test
import Polyhedra
import SCIP
import MathOptInterface as MOI

include("../src/VPolyCut.jl")
include("./helper_intersection.jl")

function run_test_simplex() 
    
    # Create a new model
    optimizer = SCIP.Optimizer()
    
    # Set Up Parameter
    #solution_path = joinpath(@__DIR__,"data","simplex.sol")
    inner = optimizer.inner
    setter = (par, val) -> SCIP.set_parameter(inner, par, val)
    turn_off_scip_miscellaneous(setter)
    turn_off_scip_heuristics(setter)
    turn_off_scip_separators(setter)
    allow_zero_power_cut(setter)
    setter("display/verblevel",5)
    
    # Add Debug Solution
    #setter("misc/debugsol",solution_path) 
    SCIP.SCIPenableDebugSol(inner.scip[]) 
 
    # The Actual Testcase
    poly1_v = Polyhedra.vrep(
        [
            [0.0,0.0],[2.5,0.0],[0,2.5]
        ]
    )
    poly1 = Polyhedra.doubledescription(poly1_v)

    # Add Variables and Constraints to the SCIP model
    n = Polyhedra.fulldim(poly1) # How many variables are there?
    x = MOI.add_variables(optimizer,n)
    for (i,x_) in enumerate(x)
        # Name the variable
        MOI.set(optimizer, MOI.VariableName(), x_,"x_$i")
    end
    hrep_to_constraint(poly1, optimizer,x)
    MOI.add_constraints(optimizer, x, MOI.Integer()) #Integrality Constraints
    
    #Set All Ones Coefficient in the Objective
    p = ones(n)
    MOI.set(
        optimizer,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.(p, x),
            0.0
        ),
    )
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)

    # Add our separator to scip
    sepa = IntersectionSeparator(scipd = inner, call_limit = 1)
    SCIP.include_sepa(inner.scip[], inner.sepas, sepa)

    # solve the problem
    SCIP.@SCIP_CALL SCIP.SCIPsolve(inner.scip[])
end

# Run The Actual Test Case
run_test_simplex()