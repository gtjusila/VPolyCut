using SCIP
using VPolyCut
using JuMP
using SCIPExperimentUtils
using Test

@testset "VPolyhedral Cut Pentagon" begin
    model = setup_scip_safe_jump_model()
    SCIPExperimentUtils.set_mode_cutting_plane_experiment(model)
    JuMP.set_attribute(model, "display/verblevel", 4)
    JuMP.set_attribute(model, "separating/gomory/freq", -1)

    @variable(model, x >= 1, Int)
    @variable(model, y >= -1, Int)
    @constraint(model, c1, x + y >= 2.6)
    @constraint(model, c2, x - y >= 1.3)
    @objective(model, Min, x) # Optimal vertex is (1.95, 0.65)

    # Include Intersection Separator
    scip = get_scip_data_from_model(model)
    VPolyCut.include_vpolyhedral_sepa(scip)

    # We use a callback to check after root node finished processing
    # Setup data dictionary
    store = Dict()
    store["scip"] = scip
    store["called"] = Ref{Bool}(false)

    # Actual Callback function 
    finishedrootnode = initiate_callback(model, store) do data
        # Mark the test case as called
        data["called"][] = true
    end

    SCIP.SCIPtransformProb(scip)
    register_callback(model, SCIP.SCIP_EVENTTYPE_NODESOLVED, finishedrootnode)
    SCIP.SCIPsolve(scip)

    # Makesure test case is called
    @test (store["called"][]) == true
end
@testset "Intersection Cut neos5" begin
    model = setup_scip_safe_jump_model()
    SCIPExperimentUtils.set_cut_selection_off(model)
    SCIPExperimentUtils.set_heuristics_emphasis_off(model)
    SCIPExperimentUtils.set_separators_emphasis_off(model)
    SCIPExperimentUtils.set_strong_branching_lookahead_off(model)
    JuMP.set_attribute(model, "display/verblevel", 4)
    JuMP.set_attribute(model, "separating/gomory/freq", -1)
    JuMP.set_attribute(model, "separating/maxroundsroot", 1)
    # Include Intersection Separator
    scip = get_scip_data_from_model(model)
    VPolyCut.include_vpolyhedral_sepa(scip)
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(scip, "neos5.mps", C_NULL)

    # We use a callback to check after root node finished processing
    # Setup data dictionary
    store = Dict()
    store["scip"] = scip
    store["called"] = Ref{Bool}(false)
    initiallpobj = -SCIP.SCIPinfinity(scip)

    firstlpsolved = initiate_callback(model, store) do data
        scip = data["scip"]
        initiallpobj = SCIP.SCIPgetLPObjval(scip)
    end

    # Actual Callback function 
    finishedrootnode = initiate_callback(model, store) do data
        scip = data["scip"]
        endofrootlpobj = SCIP.SCIPgetLPObjval(scip)
        @info initiallpobj endofrootlpobj
        # Separator should make progress towards closing LPGAP
        @test SCIP.SCIPisGT(scip, endofrootlpobj, initiallpobj) == 1
        data["called"][] = true
        SCIP.SCIPinterruptSolve(scip)
    end

    SCIP.SCIPtransformProb(scip)
    register_callback(model, SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, firstlpsolved)
    register_callback(model, SCIP.SCIP_EVENTTYPE_NODESOLVED, finishedrootnode)
    SCIP.SCIPsolve(scip)

    # Makesure test case is called
    @test (store["called"][]) == true
end