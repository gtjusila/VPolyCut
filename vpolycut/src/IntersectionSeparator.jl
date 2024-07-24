import SCIP
using JuMP
using LinearAlgebra
import MathOptInterface as MOI

DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON = false
DEBUG_PRINT_INTERSECTION_POINTS = false

"""
Intersection Cut Separator
# Fields
- scipd::SCIP.SCIPData Reference to the SCIPData object
"""
@kwdef mutable struct IntersectionSeparator <: SCIP.AbstractSeparator
    scipd::SCIP.SCIPData
    called::Int = 0
end

"""
Given a current solution, a split index, and a ray collection compute intersection points
and pararrel_rays
"""
function compute_intersection_points(current_solution, split_index, ray_collection)
    intersection_points = []
    pararrel_ray = []

    for ray in ray_collection
        low = floor(current_solution[split_index])
        up = ceil(current_solution[split_index])
        epsilon = 0

        if ray[split_index] < -EPSILON
            # Ray will hit lower bound
            epsilon = (low - current_solution[split_index]) / ray[split_index]
        elseif ray[split_index] > EPSILON
            # Ray will hit upper bound
            epsilon = (up - current_solution[split_index]) / ray[split_index]
        else
            # Ray is parallel_ray to the split
            push!(pararrel_ray, ray)
            continue
        end

        # Ray will intersect bound so add to intersection point set
        push!(intersection_points, current_solution + epsilon * ray)
    end

    return (intersection_points, pararrel_ray)
end

function modelwithsubscip()
    inner = MOI.Bridges.full_bridge_optimizer(SCIP.Optimizer(), Float64)
    model = direct_generic_model(Float64, inner)
    set_attribute(model, "display/verblevel", 0)
    backend = unsafe_backend(model)
    scip = backend.inner
    SCIP.@SCIP_CALL SCIP.SCIPsetSubscipsOff(scip, SCIP.SCIP_Bool(true))
    return model
end

function solve_separating_lp(lp_solution, intersection_points, pararrel_rays)
    dim = length(lp_solution)

    separating_lp = modelwithsubscip()

    @assert SCIP.SCIPgetSubscipsOff(unsafe_backend(separating_lp).inner) != 0

    @variable(separating_lp, -10_000 <= x[1:dim] <= 10_000)
    @variable(separating_lp, z[1:dim])

    for point in intersection_points
        new_point = point - lp_solution
        @constraint(separating_lp, sum(x[i] * new_point[i] for i = 1:dim) == 1)
    end

    for ray in pararrel_rays
        @constraint(separating_lp, sum(x[i] * ray[i] for i = 1:dim) == 0)
    end

    @constraint(separating_lp, x <= z)
    @constraint(separating_lp, -x <= z)

    @objective(separating_lp, Min, 0)

    optimize!(separating_lp)

    if is_solved_and_feasible(separating_lp)
        return Vector{SCIP.SCIP_Real}(value.(x))
    end

    return nothing
end

function project(vector::Vector{SCIP.SCIP_Real}, mask::Vector{Bool})
    return vector[mask]
end

function get_cut_vector(scip, projection_mask, projected_to_old_index, tableau, lp_sol, intersection_points, parallel_ray, dim_projected_space)
    # Project and keep projection information
    project_sol = project(lp_sol, projection_mask)
    intersection_points = [project(point, projection_mask) for point in intersection_points]
    parallel_ray = [project(ray, projection_mask) for ray in parallel_ray]

    # STEP 3: Construct and Solve Seperating LP
    separating_sol = solve_separating_lp(project_sol, intersection_points, parallel_ray)

    if isnothing(separating_sol)
        return false
    end

    b = dot(separating_sol, project_sol) + 1
    # map cut back to the original space
    full_separating_sol = zeros(get_nvars(tableau))
    for i = 1:length(separating_sol)
        full_separating_sol[projected_to_old_index[i]] = separating_sol[i]
    end
    cut_vector, b = convert_standard_row_to_general(scip, tableau, full_separating_sol, b)
    return cut_vector, b
end

"""
Construct the seperating lp return a reference to a pointer to the LPI
"""
function find_cut_from_split(
    sepa::IntersectionSeparator,
    split_index::Int64,
    lp_sol::LPSol,
    lp_rays::Vector{LPRay},
    tableau::Tableau
)::Bool
    scip = sepa.scipd
    dim = length(lp_sol)

    # STEP 3: Compute intersection Variable
    intersection_points, parallel_ray =
        compute_intersection_points(lp_sol, split_index, lp_rays)

    if DEBUG_PRINT_INTERSECTION_POINTS
        println("====================")
        println("Intersection Points")
        println(intersection_points)
        println("Parallel Ray")
        println(parallel_ray)
        println("====================")
    end

    # Project and keep projection information

    projected_to_old_index = Dict{Int,Int}()
    projection_mask = fill(false, dim)
    j = 1
    for i = 1:get_nvars(tableau)
        var = get_var_from_column(tableau, i)
        if !is_basic(var)
            projection_mask[i] = true
            projected_to_old_index[j] = i
            j += 1
        end
    end
    dim_projected_space = j - 1
    cut_vector, b = get_cut_vector(scip, projection_mask, projected_to_old_index, tableau, lp_sol, intersection_points, parallel_ray, dim_projected_space)

    separating_sol = cut_vector
    row = Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRowSepa(
        scip,
        row,
        scip.sepas[sepa],
        "",
        b,
        SCIP.SCIPinfinity(scip),
        true,
        false,
        false
    )
    vars = get_lp_variables(scip)
    infeasible = Ref{SCIP.SCIP_Bool}(0)

    for (idx, sol) in enumerate(separating_sol)
        if abs(sol) < EPSILON
            continue
        end
        SCIP.@SCIP_CALL SCIP.SCIPaddVarToRow(scip, row[], vars[idx], sol)
    end

    #SCIP.@SCIP_CALL SCIP.SCIPprintRow(scip, row[], C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPaddRow(scip, row[], true, infeasible)
    SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip, row)

    return true
end

function SCIP.exec_lp(sepa::IntersectionSeparator)
    # Aliasing for easier call
    scip = sepa.scipd
    sepa.called += 1

    # Preconditions
    @assert(SCIP.SCIPgetStage(scip) == SCIP.SCIP_STAGE_SOLVING)
    @assert(SCIP.SCIPisLPSolBasic(scip) != 0)
    @assert(SCIP.SCIPgetLPSolstat(scip) == SCIP.SCIP_LPSOLSTAT_OPTIMAL)

    # STEP 0: Get LP TAbleau data
    tableau = construct_tableau_with_constraint_matrix(scip)

    # STEP 1: Get Corner Polyhedron of current LP solution
    corner = construct_corner_polyhedron(tableau)
    lp_sol = corner.lp_sol
    lp_rays = corner.lp_rays

    if DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON
        println("====================")
        println("Seperator Called")
        println("Solution is " * string(lp_sol))
        println("LP Rays are " * string(lp_rays))
        println("====================")
    end
    # [CLEAN] Might actually not need this
    #dim = length(lp_sol)
    #if length(lp_rays) != dim
    #    return SCIP.SCIP_DIDNOTFIND
    #end

    # STEP 2: Decide Splitting Variable
    vars = get_lp_variables(scip)

    split_indices = get_all_fractional_indices(vars, 0.001)
    seperated = false

    for (i, index) in enumerate(split_indices)
        if i % 10 == 0
            println("Cut Generated $(i)")
        end
        seperated_ = find_cut_from_split(sepa, index, lp_sol, lp_rays, tableau)
        seperated = seperated || seperated_
    end

    if (seperated)
        return SCIP.SCIP_SEPARATED
    end

    return SCIP.SCIP_DIDNOTFIND
end

function include_intersection_sepa(scip::SCIP.SCIPData)
    sepa = IntersectionSeparator(scipd=scip)
    SCIP.include_sepa(scip.scip[], scip.sepas, sepa; freq=0, usessubscip=true)
end