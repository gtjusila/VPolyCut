using Test
import Polyhedra
import SCIP
import MathOptInterface as MOI

include("../src/VPolyCut.jl")
include("./utils.jl")
include("./intersection.jl")

function run_test_simplex() 
    
    # Create a new model
    optimizer = SCIP.Optimizer()
    
    # Set Up Parameter
    solution_path = joinpath(@__DIR__,"data","simplex.sol")
    inner = optimizer.inner
    setter = (par, val) -> SCIP.set_parameter(inner, par, val)
    turn_off_scip_miscellaneous(setter)
    turn_off_scip_heuristics(setter)
    turn_off_scip_separators(setter)
    allow_zero_power_cut(setter)
    setter("display/verblevel",5)
    setter("misc/debugsol",solution_path) # Add Debug Solution
 
    # The Actual Testcase
    poly1_v = Polyhedra.vrep(
        [
            [0.0,0.0],[2.5,0.0],[2.50001,2.50001],[0,2.5]
        ]
    )
    poly1 = Polyhedra.doubledescription(poly1_v)
    
    # Add Variables and Constraints to the SCIP model
    n = Polyhedra.fulldim(poly1) # How many variables are there?
    x = MOI.add_variables(optimizer,n)
    for i=1:length(x)
        MOI.set(optimizer, MOI.VariableName(), x[i],"x_$i")
    end
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
    sepa = IntersectionSeparator(scipd= inner, debug_sol_path=joinpath(@__DIR__,"simplex.sol"))
    SCIP.include_sepa(inner.scip[], inner.sepas, sepa)
    SCIP.SCIPenableVarHistory(sepa.scipd)
    # solve the problem
    SCIP.@SCIP_CALL SCIP.SCIPsolve(inner.scip[])

    @info "Number of separator calls" sepa.called
    @test sepa.called >= 1

    foo = open("transstat.txt", "w")
    file = Libc.FILE(foo) 
    SCIP.SCIPprintExpressionHandlerStatistics(sepa.scipd,file)

    # Get Vars Data
    vars = Ref{Ptr{Ptr{SCIP.SCIP_Var}}}(C_NULL)
    nvars = Ref{Cint}(0)
    SCIP.@SCIP_CALL SCIP.SCIPgetVarsData(sepa.scipd, vars, nvars, C_NULL, C_NULL, C_NULL,C_NULL)
    vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Var}},vars[],nvars[])
    foo = open("sol.txt", "w")
    file = Libc.FILE(foo) 
    SCIP.@SCIP_CALL SCIP.SCIPprintBestSol(sepa.scipd, file, 0)
    for var in vars
        bound_change = SCIP.SCIPvarGetNBdchgInfosUb(var)
        println("Bound is changed $bound_change times")
    end
end

# Run The Actual Test Case
run_test_simplex()