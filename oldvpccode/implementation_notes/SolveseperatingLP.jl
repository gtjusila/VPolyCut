function solve_separating_lp(lp_solution, intersection_points, pararrel_rays)

    dim = length(lp_solution)    
    seperating_lp = Model(SCIP.Optimizer)
    set_attribute(seperating_lp, "display/verblevel", 0)

    @variable(seperating_lp, x[1:dim])
    @variable(seperating_lp, z[1:dim])

    for point in intersection_points
        @constraint(seperating_lp, sum(x[i]*point[i] for i=1:dim) >= 1)
    end

    for ray in pararrel_rays
        @constraint(seperating_lp, sum(x[i]*ray[i] for i=1:dim) >= 0)
    end

    @constraint(seperating_lp, x  <= z)
    @constraint(seperating_lp, -x <= z)

    @objective(seperating_lp, Min, sum(z))
    
    optimize!(seperating_lp)
    #println(seperating_lp)
end