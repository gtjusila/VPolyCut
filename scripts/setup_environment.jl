using SCIP
using VPolyhedralCut.SCIPJLUtils
using Printf
using Dates
using JuMP

mutable struct ExperimentStore
    scip::SCIP.SCIPData
    result_path::String
    separator::String
    instance::String
    reference_objective::SCIP.SCIP_Real
    feasible::Bool
end

function setup_environment(execution_parameters::ExecutionParameters)
    model = setup_scip_safe_jump_model()

    set_heuristics_emphasis_off(model)
    set_separators_emphasis_off(model)
    set_cut_selection_off(model)
    set_strong_branching_lookahead_off(model)
    JuMP.set_attribute(model, "limits/restarts", 0)
    JuMP.set_attribute(model, "limits/nodes", 1)
    JuMP.set_attribute(model, "separating/maxroundsroot", 1)

    scip = get_scip_data_from_model(model)
    result_path = setup_experiment_directory(execution_parameters)
    include_separators(scip, execution_parameters.separator)
    return ExperimentStore(
        scip,
        result_path,
        execution_parameters.separator_label,
        execution_parameters.instance,
        0,
        false
    )
end

function setup_scip_object(execution_parameters::ExecutionParameters)::SCIP.SCIPData
    optimizer = SCIP.Optimizer()
    scip = optimizer.inner
    set_scip_parameters(scip, execution_parameters.easy)
    return scip
end

function setup_experiment_directory(execution_parameters::ExecutionParameters)::String
    result_path = joinpath(pwd(), "results")
    check_dir_exist_or_make(result_path)
    experiment_label = get_experiment_label(execution_parameters)
    experiment_result_path = joinpath(result_path, experiment_label)
    if (check_dir_exist_or_make(experiment_result_path))
        @error "You are already running a similar experiment"
        exit(1)
    end
    return experiment_result_path
end

function check_dir_exist_or_make(path::String)::Bool
    if !isdir(path)
        mkpath(path)
        return false
    end
    return true
end

function get_experiment_label(execution_parameters::ExecutionParameters)::String
    timestamp = get_timestamp()
    experiment_result_foldername = Printf.@sprintf(
        "%s_%s_%s",
        timestamp,
        execution_parameters.instance,
        execution_parameters.separator_label
    )
    return experiment_result_foldername
end

function get_timestamp()
    now = Dates.now()
    return Dates.format(now, "yyyymmddTHHMMSS")
end