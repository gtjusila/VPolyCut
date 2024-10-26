using SCIP
using JuMP
using VPolyhedralCut.SCIPJLUtils
using VPolyhedralCut

@kwdef mutable struct ExperimentParameters
    "Where to store the results of the experiment"
    output_path::String = ""
    "The type of separator to use"
    separator::String = ""
    "Path to the instance file"
    instance_path::String = ""
    "Is presolving on?"
    presolving::Bool = true
    "Is domain propagation on?"
    propagation::Bool = true
    "Time Limit"
    time_limit::Int = 3600
    "Zeroing Heuristic"
    zeroing_heuristic::Bool = false
end

mutable struct Experiment{T<:JuMP.AbstractModel}
    #=
    We need to redundantly keep the model and the scip pointer 
    since SCIPJLUtils work on the model level
    =#
    model::T
    scip::SCIP.SCIPData
    parameters::ExperimentParameters
    vpcsepa::Union{Nothing,VPolyhedralCut.VPCSeparator}
end

function Experiment()
    model = setup_scip_safe_jump_model()
    scip = get_scip_data_from_model(model)
    parameters = ExperimentParameters()
    return Experiment(model, scip, parameters, nothing)
end

function set_parameter(
    experiment::Experiment, parameter_name::String, value::Union{String,Integer,Bool}
)
    if parameter_name == "output_path"
        experiment.parameters.output_path = value
        return nothing
    end
    if parameter_name == "separator"
        experiment.parameters.separator = value
        return nothing
    end
    if parameter_name == "instance_path"
        experiment.parameters.instance_path = value
        return nothing
    end
    if parameter_name == "presolving" && value isa Bool
        experiment.parameters.presolving = value
        return nothing
    end
    if parameter_name == "propagation" && value isa Bool
        experiment.parameters.propagation = value
        return nothing
    end
    if parameter_name == "time_limit" && value isa Integer
        experiment.parameters.time_limit = value
        return nothing
    end
    if parameter_name == "zeroing_heuristic" && value isa Bool
        experiment.parameters.zeroing_heuristic = value
        return nothing
    end
    error("Parameter $parameter_name not found or have the wrong type")
end

function get_parameter(experiment::Experiment, parameter_name::String)
    if parameter_name == "output_path"
        return experiment.parameters.output_path
    end
    if parameter_name == "separator"
        return experiment.parameters.separator
    end
    if parameter_name == "instance_path"
        return experiment.parameters.instance_path
    end
    if parameter_name == "presolving"
        return experiment.parameters.presolving
    end
    if parameter_name == "propagation"
        return experiment.parameters.propagation
    end
    if parameter_name == "time_limit"
        return experiment.parameters.time_limit
    end
    if parameter_name == "zeroing_heuristic"
        return experiment.parameters.zeroing_heuristic
    end
    error("Parameter $parameter_name not found")
end