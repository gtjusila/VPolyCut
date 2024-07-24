# A Hack to allow calls to SCIP for julia model
function create_subscip_model()
    inner = MOI.Bridges.full_bridge_optimizer(SCIP.Optimizer(), Float64)
    model = direct_generic_model(Float64, inner)
    set_attribute(model, "display/verblevel", 0)
    backend = unsafe_backend(model)
    scip = backend.inner
    SCIP.@SCIP_CALL SCIP.SCIPsetSubscipsOff(scip, SCIP.SCIP_Bool(true))
    return model
end
