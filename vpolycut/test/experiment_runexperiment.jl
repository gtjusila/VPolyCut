using SCIP
using Dates

include("experiment_eventhandler.jl")
include("experiment_structures.jl")
include("experiment_helpers.jl")

"""
Solve `instance` with the given ExperimentConfiguration `config`.
Instance should be located on the folder test\\data\\instances
Instance should be name instance.mps and if debug is turned on solution
should be named instance.sol
"""
function run_experiment(instance::String, config::ExperimentConfiguration)
    # STEP 1: Setup logger
    main_log = open(joinpath(config.path, "experiment.log"), "w");
    println(main_log, "ExperimentType: $(config.type)")
    store = ExperimentStore()

    # STEP 2: Create a new model, Setup parameter and read problem
    optimizer = SCIP.Optimizer()
    scip::SCIP.SCIPData = optimizer.inner

    # STEP 3: Load Config to SCIP
    load_config_to_scip(scip, config)

    # STEP 4: Include V-Polyhedral Cut if needed
    if config.vpolycut
        params = IntersectionSeparatorParameter(
            call_limit = config.vpolycut_limit 
        )
        sepa = IntersectionSeparator(scipd = scip,parameter = params)
        SCIP.include_sepa(scip.scip[], scip.sepas, sepa; freq= 0, usessubscip = true)
    end

    # STEP 5: Load problem
    path = joinpath(@__DIR__,"data","instances",instance)
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(scip, path*".mps", C_NULL)

    # STEP 6: Construct and listen for event handler
    add_first_lp_event_handler(scip, main_log, store) 

    # Step 7: Load debug solution
    debug_sol = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    if config.debug
        debug_sol = load_solution(scip, debug_sol, path*".sol")
        reference_sol = SCIP.SCIPgetSolOrigObj(scip, debug_sol[])
        println(main_log,"ReferenceObjective: $(round(reference_sol,sigdigits = 4))")
        store.refobj = reference_sol
    end
  

    # STEP 8: Solve
    SCIP.@SCIP_CALL SCIP.SCIPsolve(scip)
    SCIP.@SCIP_CALL SCIP.SCIPfreeSol(scip, debug_sol)
    
    solv_stats = open(joinpath(config.path, "solving_statistic.txt"),"w")
    solv_stats_ptr = Libc.FILE(solv_stats)
    SCIP.SCIPprintStatistics(scip, solv_stats_ptr)
    
    store.finallpobj = SCIP.SCIPgetSolOrigObj(scip, C_NULL)
    store.finalgap = abs(store.finallpobj - store.refobj)

    store.gapclosed = 1 - store.finalgap/store.initialgap

    # STEP : Unload event Handlers
    #=
    SCIP.@SCIP_CALL SCIP.SCIPdropEvent(
        inner,
        SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED,
        inner.eventhdlrs[firstlp],
        C_NULL,
        C_NULL
    )
    =#
    
    println(main_log, "NFractionalVariables: $(store.nfrac)")
    println(main_log, "FirstLPObj: $(round(store.firstlpobj,sigdigits = 6))")
    println(main_log, "InitialGap: $(round(store.initialgap,sigdigits = 6))")
    println(main_log, "FinalLPObj: $(round(store.finallpobj,sigdigits = 6))")
    println(main_log, "FinalGap: $(round(store.finalgap,sigdigits = 6))")
    println(main_log, "GapClosed:  $(round(store.gapclosed,sigdigits = 6))")
    if main_log != Base.stdout
        close(main_log)
    end
end

"""
Load solution from path and place it to pointer
"""
function load_solution( scip::SCIP.SCIPData, sol::Ref{Ptr{SCIP.SCIP_Sol}}, path::String)
    partial = Ref{SCIP.SCIP_Bool}(false)
    error = Ref{SCIP.SCIP_Bool}(false)
    SCIP.@SCIP_CALL SCIP.SCIPcreateOrigSol(scip, sol, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPreadSolFile(scip, path, sol[], false, partial, error)
    if error[] == true
        @error "Error when reading solution"
    end
    if partial[] == true
        @error "Solution cannot a parial solution"
    end
    return sol
end