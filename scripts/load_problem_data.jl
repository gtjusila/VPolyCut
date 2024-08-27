function load_problem_data(experiment_store::ExperimentStore)
    load_problem_to_scip(experiment_store)
    set_reference_solution_objective(experiment_store)
end

function load_problem_to_scip(experiment_store::ExperimentStore)
    problem_path = get_problem_path(experiment_store.instance)
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(experiment_store.scip, problem_path, C_NULL)
end

function get_problem_path(instance::String)::String
    return normpath(joinpath(@__DIR__, "../instances_data/", "$(instance).mps"))
end

function set_reference_solution_objective(experiment_store::ExperimentStore)
    @assert SCIP.SCIPgetStage(experiment_store.scip) == SCIP.SCIP_STAGE_PROBLEM
    debug_sol = load_debug_solution(experiment_store)
    println(
        "Reference objective: ", SCIP.SCIPgetSolOrigObj(experiment_store.scip, debug_sol[])
    )
    experiment_store.reference_objective = SCIP.SCIPgetSolOrigObj(
        experiment_store.scip, debug_sol[]
    )
    free_debug_sol(experiment_store, debug_sol)
end

function load_debug_solution(experiment_store::ExperimentStore)::Ref{Ptr{SCIP.SCIP_Sol}}
    debug_sol = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    solution_path = get_solution_path(experiment_store.instance)
    load_solution_from_path(experiment_store.scip, debug_sol, solution_path)
    return debug_sol
end

function get_solution_path(instance::String)::String
    return normpath(joinpath(@__DIR__, "../instances_data/", "$(instance).sol"))
end

function free_debug_sol(
    experiment_store::ExperimentStore, debug_sol::Ref{Ptr{SCIP.SCIP_Sol}}
)
    SCIP.@SCIP_CALL SCIP.SCIPfreeSol(experiment_store.scip, debug_sol)
end

function load_solution_from_path(
    scip::SCIP.SCIPData, sol::Ref{Ptr{SCIP.SCIP_Sol}}, path::String
)
    partial = Ref{SCIP.SCIP_Bool}(false)
    error = Ref{SCIP.SCIP_Bool}(false)
    SCIP.@SCIP_CALL SCIP.SCIPcreateOrigSol(scip, sol, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPreadSolFile(scip, path, sol[], false, partial, error)
    if error[] == true
        @error "Error when reading solution"
    end
    if partial[] == true
        @error "Solution cannot be partial solution"
    end
end
