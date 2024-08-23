using SCIP
using Random

export initiate_callback, register_callback
mutable struct CallbackHandler <: SCIP.AbstractEventhdlr
    callback::Function
    data::Dict
end

function SCIP.eventexec(event::CallbackHandler)
    event.callback(event.data)
end

function initiate_callback(callback_function::Function, model::JuMP.AbstractModel, data::Dict)::CallbackHandler
    callback_handler = CallbackHandler(callback_function, data)
    scip = get_scip_data_from_model(model)
    # Random string is a temporary solution
    SCIP.include_event_handler(scip, callback_handler; name=randstring(10))
    return callback_handler
end

function register_callback(model::JuMP.AbstractModel, eventtype::SCIP.SCIP_EVENTTYPE, callback_handler::CallbackHandler)
    scip = get_scip_data_from_model(model)
    SCIP.catch_event(scip, eventtype, callback_handler)
end