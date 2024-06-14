import SCIP
using LinearAlgebra
using JuMP

WRITE_PATH = joinpath(pwd(),"temp")

include("../src/CornerPolyhedron.jl")
include("../src/utils.jl")
include("../src/lpi_utils.jl")

DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON = false 
DEBUG_PRINT_INTERSECTION_POINTS = false 

const LPSol = Vector{SCIP.SCIP_Real}
const LPRay = Vector{SCIP.SCIP_Real}
"""
Parameter for intersection cut separator

# Fields
- call_limit::Int64 maximum number of call to the separator (-1 for no limit). Default: -1
- file_stream::String File to write into
- epsilon::SCIP.SCIP_real Quantity to be regarded as epsilon. If negative, then epsilon will be determined at runtime from SCIP, i.e. 1000*SCIP.SCIPepsilon(). Default: -1
"""
@kwdef mutable struct IntersectionSeparatorParameter <: AbstractSeparatorParameter
    call_limit::Int64 = -1
    file_stream::IO = Base.stdout 
    epsilon::SCIP.SCIP_Real = -1.0
end

"""
Intersection Cut Separator

# Fields
- parameter::IntersectionSeparatorParameter the configuration of the separator
- called::Int64 the number of times the separator have been called
- scipd::SCIP.SCIPData Reference to the SCIPData object
"""
@kwdef mutable struct IntersectionSeparator <: SCIP.AbstractSeparator
    parameter::IntersectionSeparatorParameter = IntersectionSeparatorParameter()
    called::Int64 = 0
    scipd::SCIP.SCIPData
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
            epsilon = (low - current_solution[split_index])/ray[split_index]
        elseif ray[split_index] > EPSILON
            # Ray will hit upper bound
            epsilon = (up - current_solution[split_index])/ray[split_index]
        else
            # Ray is parallel_ray to the split
            push!(pararrel_ray,ray)
            continue
        end

        # Ray will intersect bound so add to intersection point set
        push!(intersection_points, current_solution + epsilon*ray)
    end

    return (intersection_points,pararrel_ray)
end

"""
Construct the seperating lp return a reference to a pointer to the LPI
"""
function constructSeperatingLP(lp_sol, intersection_points, parallel_ray)
    dim = length(lp_sol)
    lpi = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(lpi, C_NULL, "Seperating LP", SCIP.SCIP_OBJSEN_MINIMIZE)

    # Add Variables
    SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(lpi[],dim,zeros(dim),-10000*ones(dim),10000*ones(dim),C_NULL, 0,C_NULL, C_NULL, C_NULL )
    
    # Point constraint
    for (idx,point) in enumerate(intersection_points)
        @assert(length(point) == dim)
        add_row_to_lpi(lpi[]; alpha = (point-lp_sol), lhs = 1.0, rhs = 1.0, row_name = "point"*string(idx))
    end

    # Ray constraint
    for (idx,ray) in enumerate(parallel_ray)
        @assert(length(ray) == dim)
        add_row_to_lpi(lpi[]; alpha = ray, lhs = EPSILON, row_name = "ray"*string(idx))
    end

    SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(lpi[],dim,ones(dim),-10000*ones(dim),10000*ones(dim),C_NULL, 0,C_NULL, C_NULL, C_NULL)
    
    for i=1:dim
        alpha = zeros(SCIP.SCIP_Real, 2*dim)

        alpha[i] = -1.0
        alpha[i+dim] = -1.0
        add_row_to_lpi(lpi[]; alpha = alpha, rhs = 0.0)

        alpha[i] = 1.0
        alpha[i+dim] = -1.0
        add_row_to_lpi(lpi[]; alpha = alpha, rhs = 0.0) 
    end

    SCIP.@SCIP_CALL SCIP.SCIPlpiWriteLP(lpi[],"test.lp")
    return lpi
end

function find_cut_from_split(sepa::IntersectionSeparator, split_index::Int64, lp_sol::LPSol, lp_rays::Vector{LPRay})::Bool
    scip = sepa.scipd
    dim = length(lp_sol)
    # STEP 3: Compute intersection Variable
    intersection_points, parallel_ray = compute_intersection_points(lp_sol, split_index, lp_rays)
    
    if DEBUG_PRINT_INTERSECTION_POINTS
        println("====================")
        println("Intersection Points")
        println(intersection_points)
        println("Parallel Ray")
        println(parallel_ray)
        println("====================")
    end

    # STEP 3: Construct and Solve Seperating LP
    lpi = constructSeperatingLP(lp_sol, intersection_points, parallel_ray) 

    try
        solve_lpi(lpi[])
    catch
        return false 
    end

    separating_sol = get_lpi_solution_vector(lpi[])
    SCIP.@SCIP_CALL SCIP.SCIPlpiFree(lpi)
    
    separating_sol = separating_sol[1:dim]
    b = dot(separating_sol,lp_sol) + 1

    row = Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRowSepa(scip, row, scip.sepas[sepa], "", b,SCIP.SCIPinfinity(scip), true, false, false)      
    vars = get_lp_variables(scip)
    infeasible = Ref{SCIP.SCIP_Bool}(0)

    for (idx, sol) in enumerate(separating_sol)
        if abs(sol) < EPSILON continue end
        SCIP.@SCIP_CALL SCIP.SCIPaddVarToRow(scip, row[], vars[idx],sol)
    end
    #SCIP.@SCIP_CALL SCIP.SCIPprintRow(scip, row[], C_NULL)

    SCIP.@SCIP_CALL SCIP.SCIPaddRow(scip,row[],true, infeasible)  
    SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip,row)
    
    return true
end
# wrapper for writing to log files
function SCIP.exec_lp(sepa::IntersectionSeparator)
    redirect_stdout(sepa.parameter.file_stream) do 
        _exec_lp(sepa)
    end
end

function _exec_lp(sepa::IntersectionSeparator)
    # Aliasing for easier call
    scip = sepa.scipd

    # Preconditions
    @assert(SCIP.SCIPgetStage(scip) == SCIP.SCIP_STAGE_SOLVING)
    @assert(SCIP.SCIPisLPSolBasic(scip) != 0)
    @assert(SCIP.SCIPgetLPSolstat(scip) == SCIP.SCIP_LPSOLSTAT_OPTIMAL) 

    # Only limit 1 call 
    if sepa.parameter.call_limit > 0  && sepa.called >= sepa.parameter.call_limit
        return SCIP.SCIP_DIDNOTRUN
    end
    sepa.called +=1 
 
    # STEP 1: Get Corner Polyhedron of current LP solution
    lp_sol = nothing
    lp_rays = nothing
    lp_sol, lp_rays = get_corner_polyhedron(scip)
    
    if DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON
        println("====================")
        println("Seperator Called")
        println("Solution is "*string(lp_sol[abs.(lp_sol) .> EPSILON]))
        println("LP Rays are "*string(lp_rays))
        println("====================")
    end

    # [CLEAN] Might actually not need this
    dim = length(lp_sol)
    if length(lp_rays) != dim
        return SCIP.SCIP_DIDNOTFIND 
    end
    
    # STEP 2: Decide Splitting Variable
    vars = get_lp_variables(scip)
    split_index = get_most_fractional_index(vars)
    
    if split_index == -1 
        @warn "No Splitting Variable"
        return SCIP.SCIP_DIDNOTFIND
    end
    
    println(lp_sol[split_index])
    split_indices = get_all_fractional_indices(vars,0.001)
    
    @info length(split_indices)
    seperated = false

    for (i,index) in enumerate(split_indices)
        if i%10 == 0 println("Cut Generated $(i)") end
        seperated_ = find_cut_from_split(sepa, index, lp_sol, lp_rays)
        seperated = seperated || seperated_
    end

    if (seperated)
        return SCIP.SCIP_SEPARATED
    end

    return SCIP.SCIP_DIDNOTFIND    
end
