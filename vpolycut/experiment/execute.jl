#
# execute.jl
# The code containing a single run of the experiment
#
import SCIP

include("helpers.jl")
include("../src/IntersectionSeperator.jl")

function execute(settings::Dict)
    # STEP 0: Setup
    optimizer = SCIP.Optimizer()
    scip::SCIP.SCIPData = optimizer.inner
    setter_param = (par, val) -> SCIP.set_parameter(scip, par, val)
    setscipsettings(setter_param)

    # STEP 1: Activate desired cut
    settings["mode"] = strip(settings["mode"])
    if settings["mode"] == "gomory"
        includegomorysepa(setter_param)
    elseif settings["mode"] == "vpc"
        includevpcsepa(scip)
    else
        println("Unrecognized Experiment Mode. Terminating.")
        exit(1)
    end

    # STEP 2: Load problem
    instance = strip(settings["instance"])
    path = joinpath(abspath("./data"), instance)
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(scip, path * ".mps", C_NULL)

    # STEP 3: Load Debugging solution
    debug_sol = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    load_solution(scip, debug_sol, path * ".sol")
    reference_sol = SCIP.SCIPgetSolOrigObj(scip, debug_sol[])
    println("ReferenceObjective: $(round(reference_sol,sigdigits = 4))")
    SCIP.@SCIP_CALL SCIP.SCIPfreeSol(scip, debug_sol)

    # STEP 4: Solve
    SCIP.@SCIP_CALL SCIP.SCIPsolve(scip)

    initialdualbound = SCIP.SCIPgetFirstLPDualboundRoot(scip)
    println("Initial Dual Bound $(initialdualbound)")
    finaldualbound = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    println("Final Dual Bound $(finaldualbound)")

    initialgap = abs(initialdualbound - reference_sol)
    finalgap = abs(finaldualbound - reference_sol)

    gapclosed = 1 - finalgap / initialgap

    solv_stats = open("solving_statistic.txt", "w")
    solv_stats_ptr = Libc.FILE(solv_stats)
    SCIP.SCIPprintStatistics(scip, solv_stats_ptr)

    println("Gap Closed $(gapclosed)")
end