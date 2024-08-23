# create a scip safe jump model
using JuMP

export setup_scip_safe_jump_model, get_scip_data_from_model

function setup_scip_safe_jump_model()
    inner = MOI.Bridges.full_bridge_optimizer(SCIP.Optimizer(), Float64)
    model = JuMP.direct_generic_model(Float64, inner)
    setter = get_parameter_setter_function(model)

    # For stability
    setter("limits/restarts", 0)

    return model
end

function get_scip_data_from_model(model::T) where {T<:JuMP.AbstractModel}
    backend = JuMP.unsafe_backend(model)
    return backend.inner
end
