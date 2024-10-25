using SCIP
using VPolyhedralCut.SCIPJLUtils
using VPolyhedralCut
using JuMP

function load_problem_to_scip(experiment::Experiment)
    instance_relative_path = get_parameter(experiment, "instance_path")
    instance_path = abspath(instance_relative_path)
    @info "Loading instance from $instance_path"
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(experiment.scip, instance_path, C_NULL)
end

function setup_scip_parameter(experiment::Experiment)
    # The following settings are always applied
    set_heuristics_emphasis_off(experiment.model)
    set_separators_emphasis_off(experiment.model)
    set_cut_selection_off(experiment.model)
    set_strong_branching_lookahead_off(experiment.model)

    JuMP.set_attribute(experiment.model, "limits/restarts", 0)
    JuMP.set_attribute(experiment.model, "limits/nodes", 1)
    JuMP.set_attribute(
        experiment.model, "limits/time", get_parameter(experiment, "time_limit")
    )
    JuMP.set_attribute(experiment.model, "separating/maxroundsroot", 1)

    # The following settings are applied optionaly 
    if get_parameter(experiment, "presolving") == false
        set_presolving_emphasis_off(experiment.model)
    end
    if get_parameter(experiment, "propagation") == false
        set_root_node_propagation_off(experiment.model)
    end

    # Include separator
    include_separator(experiment)
end

function include_separator(experiment::Experiment)
    separator = get_parameter(experiment, "separator")
    if separator == "gomory"
        include_gomory_separator(experiment)
        return nothing
    end
    if separator == "intersection"
        VPolyhedralCut.include_intersection_sepa(experiment.scip)
        return nothing
    end
    if separator == "vpc"
        VPolyhedralCut.include_vpolyhedral_sepa(experiment.scip; n_leaves=64)
        return nothing
    end
    error("Separator $separator not found")
end

function include_gomory_separator(experiment::Experiment)
    setter = (par, val) -> SCIP.set_parameter(experiment.scip, par, val)
    setter("separating/gmi/freq", 0)
    setter("separating/gmi/priority", 9999)
    setter("separating/gmi/maxsuppabs", 5000)
    setter("separating/gmi/dynamiccuts", false)
    setter("separating/gmi/maxsupprel", 1.0)
    setter("separating/gmi/forcecuts", true)
end

#=
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
=#