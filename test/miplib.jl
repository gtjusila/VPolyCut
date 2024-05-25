using Test
import Polyhedra
import SCIP
import MathOptInterface as MOI

include("../utils.jl")
include("../intersection.jl")

@testset "MIPLIB" begin 

    # STEP 1: Create a new model, Setup parameter and read problem
    optimizer = SCIP.Optimizer()
    inner = optimizer.inner
    sepa_set_scip_parameters_with_presolve((par, val) -> SCIP.set_parameter(inner, par, val))
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(inner.scip[],joinpath(@__DIR__,"mps","gen-ip054.mps"),C_NULL)
    SCIP.set_parameter(inner,"misc/debugsol",joinpath(@__DIR__,"mps","gen-ip054.sol"))
    SCIP.SCIPenableDebugSol(inner.scip[]) 
    # Add our separator to scip
    sepa = IntersectionSeparator(scipd = inner, debug_sol_path = joinpath(@__DIR__,"mps","gen-ip054.sol"))
    SCIP.include_sepa(inner.scip[], inner.sepas, sepa; freq= 0) 
    
    # Do the actual solve
    SCIP.@SCIP_CALL SCIP.SCIPsolve(inner.scip[])

    #Write stats
    foo = open("stat.txt", "w")
    file = Libc.FILE(foo) 
    #SCIP.SCIPprintSeparatorStatistics(sepa.scipd,file)
    
    vars = Ref{Ptr{Ptr{SCIP.SCIP_VAR}}}(C_NULL)
    nvars = Ref{Cint}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPgetVarsData(inner.scip[],vars,nvars,C_NULL,C_NULL,C_NULL,C_NULL)
    vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Var}},vars[],nvars[])
    nvars = nvars[]
    solution = zeros(nvars)
    for i=1:nvars
        solution[i] = SCIP.SCIPvarGetSol(vars[i], 0)
    end
    #println(solution)
end