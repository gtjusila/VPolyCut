using SCIP
using Dates
using Printf

include("helper_parametersetter.jl")
include("helper_eventhandler.jl")

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
- debug::Bool use debug solution? Default: False
- verbosity::Int Display verbosity level. Default: 5
- log_file::string Directory to write the log file. Default: Base.stdout 
- vpolycut::Bool should V Polyhedral Cut be used? Default: True
- vpolycut_limit::Int64 Maximum number of times the V Polyhedral Cut Seperator can be called? (-1 for no limit) Default: -1
- gomory::Bool Override seperator settings and turn on Gomory Cut for the root node? Default: False 
- node_limit::Int SCIP node limit. (-1 for no limit) Default: -1
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
    vpolycut_limit::Int64 = -1
    gomory::Bool = false
    node_limit::Int64 = -1
end

"""
Solve `instance` with the given ExperimentConfiguration `config`.
Instance should be located on the folder test\\data\\instances
Instance should be name instance.mps and if debug is turned on solution
should be named instance.sol
"""
function run_experiment(instance::String, config::ExperimentConfiguration)
    # STEP 1: Configure experiment folder
    experiment_path = setup_experiment_folder(instance)

    # STEP 2: Setup logger
    main_log = nothing;
    try
        main_log = open(joinpath(experiment_path, "experiment.log"), "w+");
    catch
        main_log = Base.stdout
    end

    # STEP 3: Create a new model, Setup parameter and read problem
    optimizer = SCIP.Optimizer()
    inner::SCIP.SCIPData = optimizer.inner

    # STEP 4: Set parameter
    setter = (par, val) -> SCIP.set_parameter(inner, par, val)
    if !config.separator turn_off_scip_separators(setter) end
    if !config.heuristics turn_off_scip_heuristics(setter) end
    if !config.presolve setter("presolving/maxrounds", 0) end
    if !config.propagation setter("propagating/maxroundsroot",0) end
    if !config.conflict setter("conflict/enable",false) end
    if config.zero_cut allow_zero_power_cut(setter) end
    if !config.symmetry setter("misc/usesymmetry", 0) end
    if !config.gomory setter("separating/gomory/freq", 0) end
    setter("display/verblevel",config.verbosity) 
    setter("limits/nodes",config.node_limit)

    # STEP 5: Load problem with debug solution (if debug is turned on)
    path = joinpath(@__DIR__,"data","instances",instance)
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(inner, path*".mps", C_NULL)
    
    # Commented this part out we do manual check instead
    # if config.debug
    #    SCIP.set_parameter(inner,"misc/debugsol",joinpath(@__DIR__,"data","instances",instance*".sol"))
    #    SCIP.SCIPenableDebugSol(inner)
    # end

    # STEP 6: Construct and listen for event handler
    # SCIP must be in problem stage when we include event handler 
    # and in transformed stage when we call catch event
    firstlp = RecordSolutionAfterFirstRootLPSolve(scipd = inner)
    SCIP.include_event_handler(inner.scip[], inner.eventhdlrs, firstlp)
    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(inner)
    SCIP.@SCIP_CALL SCIP.SCIPcatchEvent(
        inner.scip[], 
        SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, 
        inner.eventhdlrs[firstlp],
        C_NULL,
        C_NULL
    )

    # Step 7: Load debug solution
    debug_sol = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    if config.debug
        partial = Ref{SCIP.SCIP_Bool}(false)
        error = Ref{SCIP.SCIP_Bool}(false)
        SCIP.@SCIP_CALL SCIP.SCIPcreateOrigSol(inner,debug_sol, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPreadSolFile(inner, path*".sol", debug_sol[], false, partial, error)
        if error[] == true 
            @error "Error when reading solution"
        end
        if partial[] == true
            @error "Solution cannot a parial solution"
        end
        reference_sol = SCIP.SCIPgetSolOrigObj(inner, debug_sol[])
        println(main_log,"Reference Objective is $(round(reference_sol,sigdigits = 4))")
        flush(main_log)
    end

    
    # STEP 8: Include V-Polyhedral Cut
    if config.vpolycut
        params = IntersectionSeparatorParameter(
            call_limit = vpolycut_limit 
        )
        sepa = IntersectionSeparator(scipd = inner,parameter = params)
        SCIP.include_sepa(inner.scip[], inner.sepas, sepa; freq= 0, usessubscip = true)
    end

    # STEP 9: Solve
    SCIP.@SCIP_CALL SCIP.SCIPsolve(inner)

    # STEP 10: Unload event Handlers
    #=
    SCIP.@SCIP_CALL SCIP.SCIPdropEvent(
        inner,
        SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED,
        inner.eventhdlrs[firstlp],
        C_NULL,
        C_NULL
    )
    =#
    if main_log != Base.stdout
        close(main_log)
    end
end

"""
Setup experiment folder. The format is xxxxx_instancename_date where xxxxx is a random
alphanumeric code
"""
function setup_experiment_folder(instance::String)
    if !isdir("experiments")
        mkdir("experiments")
    end
    code = lowercase(randstring(5))
    date = Dates.format(Dates.today(), "yyyymmdd")
    path = joinpath(pwd(), "experiments", "$(code)_$(instance)_$(date)")
    mkdir(path)
    println("Created Experiment Folder $(path)")
    return path
end
