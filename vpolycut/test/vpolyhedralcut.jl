using SCIP
using VPolyCut
using JuMP
using SCIPExperimentUtils
using Test

@testset "VPolyhedral Cut Pentagon" begin
    model = setup_scip_safe_jump_model()
    SCIPExperimentUtils.set_mode_cutting_plane_experiment(model)
    JuMP.set_attribute(model, "display/verblevel", 4)

    @variable(model, 0 <= x <= 2, Int)
    @variable(model, 0 <= y <= 2, Int)
    @constraint(model, c1, x + y <= 2.5)
    @constraint(model, c2, x - y <= 1.3)
    @objective(model, Min, -x) # Optimal vertex is (1.9, 0.6)

    # Include Intersection Separator
    scip = get_scip_data_from_model(model)
    VPolyCut.include_vpolyhedral_sepa(scip)

    # Originally solution should be feasible
    vars = SCIP.SCIPgetVars(scip)
    nvars = SCIP.SCIPgetNVars(scip)
    vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, vars, nvars)
    sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 1.9)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 0.6)
    feasible = Ref{SCIP.SCIP_Bool}(0)

    SCIP.@SCIP_CALL SCIP.SCIPcheckSol(
        scip, sol[], 0, 1, 1, 0, 1, feasible
    )
    @test feasible[] == 1
    SCIP.@SCIP_CALL SCIP.SCIPfreeSol(scip, sol)

    # We use a callback to check after root node finished processing
    # Setup data dictionary
    store = Dict()
    store["scip"] = scip
    store["called"] = Ref{Bool}(false)

    # Actual Callback function 
    finishedrootnode = initiate_callback(model, store) do data
        scip = data["scip"]

        # Catch when root node is branched or solved
        # Initial LP solution should have been cut off by now 
        feasible = Ref{SCIP.SCIP_Bool}(0)
        sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 1.9)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 0.6)
        SCIP.@SCIP_CALL SCIP.SCIPcheckSol(
            scip, sol[], false, 1, 1, 0, 1, feasible
        )
        @test feasible[] == 0
        SCIP.@SCIP_CALL SCIP.SCIPfreeSol(scip, sol)
        SCIP.SCIPinterruptSolve(scip)

        # Mark the test case as called
        data["called"][] = true
    end

    SCIP.SCIPtransformProb(scip)
    register_callback(model, SCIP.SCIP_EVENTTYPE_NODESOLVED, finishedrootnode)
    SCIP.SCIPsolve(scip)

    # Makesure test case is called
    @test (store["called"][]) == true
end
