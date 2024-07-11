#
# execute.jl
# The code containing a single run of the experiment
#
import SCIP

include("setscipparam.jl")
include("helpers.jl")
include("writeoutput.jl")
include("../src/IntersectionSeparator.jl")

function execute(settings::Dict)
    # STEP 0: Setup
    optimizer = SCIP.Optimizer()
    scip::SCIP.SCIPData = optimizer.inner
    setter_param = (par, val) -> SCIP.set_parameter(scip, par, val)
    setscipsettings(setter_param, settings["easy"])
    settings["mode"] = string(strip(settings["mode"]))
    settings["instance"] = string(strip(settings["instance"]))
    result_path = setupexperimentdirectory(settings["instance"], settings["mode"])

    # STEP 1: Activate desired cut
    if settings["mode"] == "gomory"
        includegomorysepa(setter_param)
    elseif settings["mode"] == "vpc"
        includeintersectionsepa(scip)
    else
        println("Unrecognized Experiment Mode. Terminating.")
        exit(1)
    end

    # STEP 2: Load problem
    path = joinpath(abspath("./data"), settings["instance"])
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(scip, path * ".mps", C_NULL)

    # STEP 3: Load Debugging solution for reference_obj
    debug_sol = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    load_solution(scip, debug_sol, path * ".sol")
    reference_obj = SCIP.SCIPgetSolOrigObj(scip, debug_sol[])
    SCIP.@SCIP_CALL SCIP.SCIPfreeSol(scip, debug_sol)

    # STEP 4: Solve
    SCIP.@SCIP_CALL SCIP.SCIPsolve(scip)

    # STEP 5: Check if debug solution is still feasible:
    debug_sol = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    feasible = Ref{SCIP.SCIP_Bool}(C_NULL)
    load_solution(scip, debug_sol, path * ".sol")
    SCIP.@SCIP_CALL SCIP.SCIPcheckSol(scip, debug_sol[], true, true, true, true, true, feasible)
    feasible = feasible[]

    # STEP 5: Check if debug solution is still feasible:
    debug_sol = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    feasible = Ref{SCIP.SCIP_Bool}(C_NULL)
    load_solution(scip, debug_sol, path * ".sol")
    SCIP.@SCIP_CALL SCIP.SCIPcheckSol(scip, debug_sol[], true, true, true, true, true, feasible)
    feasible = feasible[]


    writeoutput(result_path, scip, settings, reference_obj, feasible)
end