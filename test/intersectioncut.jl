using SCIP
using VPolyhedralCut
using VPolyhedralCut.SCIPJLUtils
using JuMP
using Test

# Test the Intersection Cut over a Simplex
@testset "Intersection Simplex" begin
    model = setup_scip_safe_jump_model()
    set_mode_cutting_plane_experiment(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    @variable(model, x >= 0, Int)
    @variable(model, y >= 0, Int)
    @constraint(model, c1, x + y <= 1.5)
    @objective(model, Min, -x)

    # Include Intersection Separator
    scip = get_scip_data_from_model(model)
    VPolyhedralCut.include_intersection_sepa(scip)

    # Originally solution should be feasible
    vars = SCIP.SCIPgetVars(scip)
    nvars = SCIP.SCIPgetNVars(scip)
    vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, vars, nvars)
    sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 1.5)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 0.0)
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
        feasible = Ref{SCIP.SCIP_Bool}(0)
        sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 1.5)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 0.0)
        SCIP.@SCIP_CALL SCIP.SCIPcheckSol(
            scip, sol[], true, 1, 1, 0, 1, feasible
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

# Test the Intersection Cut over a Box constraint
@testset "Intersection Box" begin
    model = setup_scip_safe_jump_model()
    set_mode_cutting_plane_experiment(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    @variable(model, x >= 0, Int)
    @variable(model, y >= 0, Int)
    @constraint(model, c1, x + 10e-8 * y <= 1.6)
    @constraint(model, c2, 10e-8 * x + y <= 1.6)
    @objective(model, Min, -x)

    scip = get_scip_data_from_model(model)
    VPolyhedralCut.include_intersection_sepa(scip)

    # Makesure the original LP Solution is cut off 

    # Originally solution should be feasible
    vars = SCIP.SCIPgetVars(scip)
    nvars = SCIP.SCIPgetNVars(scip)
    vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, vars, nvars)
    sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 1.58)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 1.58)
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
        # [1.5 1] should be LP infeasible here
        feasible = Ref{SCIP.SCIP_Bool}(0)
        sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 1.58)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 1.58)
        SCIP.@SCIP_CALL SCIP.SCIPcheckSol(
            scip, sol[], true, 1, 1, 0, 1, feasible
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

# Test the Intersection Cut over a Box constraint
@testset "Intersection Reverse" begin
    println("Testcase reverse")
    model = setup_scip_safe_jump_model()
    set_mode_cutting_plane_experiment(model)
    JuMP.set_attribute(model, "display/verblevel", 0)

    @variable(model, x <= 2.5, Int)
    @variable(model, 0 <= y, Int)
    @constraint(model, c1, x + 10e-8 * y >= 0.5)
    @constraint(model, c2, 10e-8 * x + y <= 2.5)
    @objective(model, Min, x)

    scip = get_scip_data_from_model(model)
    VPolyhedralCut.include_intersection_sepa(scip)

    # Makesure the original LP Solution is cut off 

    # Originally solution should be feasible
    vars = SCIP.SCIPgetVars(scip)
    nvars = SCIP.SCIPgetNVars(scip)
    vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, vars, nvars)
    sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 0.55)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 2.45)
    feasible = Ref{SCIP.SCIP_Bool}(0)
    # [1.5 1] should be LP feasible here
    SCIP.@SCIP_CALL SCIP.SCIPcheckSol(
        scip, sol[], 1, 1, 1, 0, 1, feasible
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
        # [1.5 1] should be LP infeasible here
        feasible = Ref{SCIP.SCIP_Bool}(0)
        sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 0.55)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 2.45)
        SCIP.@SCIP_CALL SCIP.SCIPcheckSol(
            scip, sol[], true, 1, 1, 0, 1, feasible
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

@testset "Intersection Cut pentagon" begin
    model = setup_scip_safe_jump_model()
    set_mode_cutting_plane_experiment(model)
    JuMP.set_attribute(model, "display/verblevel", 4)

    @variable(model, 0 <= x <= 2, Int)
    @variable(model, 0 <= y <= 2, Int)
    @constraint(model, c1, x + y <= 2.5)
    @constraint(model, c2, x - y <= 1.3)
    @objective(model, Min, -x)

    # Include Intersection Separator
    scip = get_scip_data_from_model(model)
    VPolyhedralCut.include_intersection_sepa(scip)

    # Originally solution should be feasible
    vars = SCIP.SCIPgetVars(scip)
    nvars = SCIP.SCIPgetNVars(scip)
    vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, vars, nvars)
    sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 1)
    SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 0.0)
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
        feasible = Ref{SCIP.SCIP_Bool}(0)
        sol = Ref{Ptr{SCIP.SCIP_SOL}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPcreateSol(scip, sol, C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[1], 1.5)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(scip, sol[], vars[2], 1.6)
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

@testset "Intersection Cut neos5" begin
    model = setup_scip_safe_jump_model()
    set_cut_selection_off(model)
    set_heuristics_emphasis_off(model)
    set_separators_emphasis_off(model)
    set_strong_branching_lookahead_off(model)
    JuMP.set_attribute(model, "display/verblevel", 4)
    JuMP.set_attribute(model, "separating/gomory/freq", -1)
    JuMP.set_attribute(model, "separating/maxroundsroot", 1)
    # Include Intersection Separator
    scip = get_scip_data_from_model(model)
    VPolyhedralCut.include_intersection_sepa(scip)
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(scip, "../instances_data/neos5.mps", C_NULL)

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
