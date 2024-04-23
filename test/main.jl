using Test
import Polyhedra
import SCIP
import MathOptInterface as MOI

include("../src/VPolyCut.jl")
include("./utils.jl")

# Setup Temp Forlders For Experiment
rm(joinpath("temp"), force=true, recursive=true)
mkdir(joinpath(pwd(),"temp"))

@testset "Simplex" begin
    # Create a new model
    optimizer = SCIP.Optimizer()
    inner = optimizer.inner
    sepa_set_scip_parameters((par, val) -> SCIP.set_parameter(inner, par, val))    

    # The Actual Testcase
    poly1_v = Polyhedra.vrep([[0.0,0.0],[2.25,0],[0,2.25]])
    poly1 = Polyhedra.doubledescription(poly1_v)
    
    # Add Variables and Constraints to the SCIP model
    n = Polyhedra.fulldim(poly1) # How many variables are there?
    x = MOI.add_variables(optimizer, n)
    hrep_to_constraint(poly1, optimizer,x)
    MOI.add_constraints(optimizer, x, MOI.Integer())

    #Add All Ones Coefficient in the Objective
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
    sepa = VPolySeparator(inner)
    SCIP.include_sepa(inner.scip[], inner.sepas, sepa)

    # solve the problem
    SCIP.@SCIP_CALL SCIP.SCIPsolve(inner.scip[])

    @info "Number of separator calls" sepa.called
    @test sepa.called >= 1
end