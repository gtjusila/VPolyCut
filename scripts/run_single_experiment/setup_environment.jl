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
        experiment.vpcsepa = VPolyhedralCut.include_vpolyhedral_sepa(
            experiment.scip;
            n_leaves=get_parameter(experiment, "number_of_leaves"),
            write_log=!get_parameter(experiment, "use_stdout"),
            log_directory=get_parameter(experiment, "output_path"),
            lp_solving_method=get_parameter(experiment, "lp_solving_method")
        )
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