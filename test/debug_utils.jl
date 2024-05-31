#[CLEAN] Debug Function

function check_in_corner_polyhedron(vertex, rays, point)
    model = Model(SCIP.Optimizer)
    set_attribute(model, "display/verblevel",0)
    @variable(model, l[i=1:length(rays)]>=-0.1)
    @variable(model, z)
    @constraint(model, l .>= z)
    @objective(model, Max, z)
    R = hcat(rays...)
    @constraint(model, R*l + vertex .== point)
    optimize!(model)
    if !is_solved_and_feasible(model)
        @warn "Reference Solution Not In Corner Polyhedron"
    else
        println("Solution Is In The Corner Polyhedron")
    end
end

function verify_reference_solution_in_disjunction(sepa::IntersectionSeparator,reference_solution, split_index)
    lp_cols = get_lp_columns(sepa.scipd)
    col_num = length(lp_cols)
    split_var = SCIP.SCIPcolGetVar(lp_cols[split_index])
    split_var_value = SCIP.SCIPvarGetSol(split_var,1)
    reference = SCIP.SCIPgetSolVal(sepa.scipd, reference_solution[], split_var)
    println(string(split_var_value)*" vs "*string(reference))
    
    found = false
    SCIP.@SCIP_CALL SCIP.SCIPstartProbing(sepa.scipd)
    SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(sepa.scipd)
    SCIP.@SCIP_CALL SCIP.SCIPchgVarUbProbing(sepa.scipd, split_var, floor(split_var_value))
    println("Testing For Lower Than")
    feasible = Ref{SCIP.SCIP_Bool}(0)
    SCIP.@SCIP_CALL SCIP.SCIPcheckSol(sepa.scipd,reference_solution[], 1,1,1,0,1,feasible)
    println(feasible[]==1)
    found = found || (feasible[]==1)
    SCIP.@SCIP_CALL SCIP.SCIPbacktrackProbing(sepa.scipd, 0)
    SCIP.@SCIP_CALL SCIP.SCIPnewProbingNode(sepa.scipd)
    SCIP.@SCIP_CALL SCIP.SCIPchgVarLbProbing(sepa.scipd,split_var, ceil(split_var_value))
    println("Testing For Larger Than")
    feasible = Ref{SCIP.SCIP_Bool}(0)
    SCIP.@SCIP_CALL SCIP.SCIPcheckSol(sepa.scipd,reference_solution[], 1,1,1,0,1,feasible)
    println(feasible[]==1)
    found = found || (feasible[]==1) 
    SCIP.@SCIP_CALL SCIP.SCIPendProbing(sepa.scipd)
end
