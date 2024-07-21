# 
# This file contains the code to retrieve the corner polyhedron from the current SCIP Pointer 
# along with the helper function to get Tableau Information
#
import SCIP

const LPRay = Vector{SCIP.SCIP_Real}
const LPSol = Vector{SCIP.SCIP_Real}
struct CornerPolyhedron
    lp_sol::LPSol
    lp_rays::Vector{LPRay}
end

"""
Returns the pair (lp_sol, lp_rays) where lp_sol is the current solution and lp_rays is a vector
containing lp_rays each is of the form Vector{SCIP.SCIP_Real}
"""
function get_corner_polyhedron(tableau::LPTableau)::CornerPolyhedron
    # Initiate a vector to collect corner polyhedron ray
    dim = get_num_tableau_objects(tableau)
    rays = get_non_basic_rays(tableau)
    sol = get_solution_vector(tableau)

    projection_mask = fill(false, dim)
    for i = 1:get_num_lp_cols(tableau)
        projection_mask[i] = true
    end

    # Project only to the space contained by LP Columns 
    sol = project(sol, projection_mask)
    rays = [project(ray, projection_mask) for ray in rays]

    return CornerPolyhedron(sol, rays)
end

function project(vector::Vector{SCIP.SCIP_Real}, mask::Vector{Bool})
    return vector[mask]
end

function get_solution_vector(tableau::LPTableau)::Vector{SCIP.SCIP_Real}
    dim = get_num_tableau_objects(tableau)
    solution = zeros(dim)

    for i = 1:dim
        col = get_column_object(tableau, i)
        solution[i] = get_solution(col)
    end

    return solution
end

function get_solution(col::LPColumn)::SCIP.SCIP_Real
    return get_col_sol(col)
end

function get_solution(row::LPRow)::SCIP.SCIP_Real
    return -get_row_slack(row)
end

function get_non_basic_rays(tableau::LPTableau)::Vector{LPRay}
    ray_collection = Vector{LPRay}(undef, 0)

    for i = 1:get_num_tableau_objects(tableau)
        col = get_column_object(tableau, i)
        if get_basis_status(col) == SCIP.SCIP_BASESTAT_BASIC
            continue
        end
        ray = construct_non_basic_ray(tableau, col)
        push!(ray_collection, ray)
    end

    return ray_collection
end

function construct_non_basic_ray(tableau::LPTableau, col::LPObject)::LPRay
    direction = 1.0

    if (get_basis_status(col) == SCIP.SCIP_BASESTAT_UPPER)
        direction = -1.0
    elseif (get_basis_status(col) == SCIP.SCIP_BASESTAT_LOWER)
        direction = 1.0
    elseif (get_basis_status(col) == SCIP.SCIP_BASESTAT_ZERO)
        #Safekeeping: Should Never Happen unless polyhedron is not pointed
        error("SCIP_BASESTAT_ZERO encountered")
    else
        error("SCIP_BASESTAT Status Undefined")
    end

    # SCIP assumes that the constraint matrix is in the form [A I] where I
    # are columns corresponding to the slack variables. Hence, if the 
    # nonbasic column is a slack variable, then the direction is reversed
    if isa(col, LPRow)
        direction = -direction
    end

    # Construct ray r
    dim = get_num_tableau_objects(tableau)
    ray = zeros(dim)

    col_idx = get_object_column(tableau, col)
    ray[col_idx] = direction

    for row_idx = 1:get_num_basic_objects(tableau)
        row_obj = get_row_object(tableau, row_idx)
        row_obj_col = get_object_column(tableau, row_obj)
        ray[row_obj_col] = -direction * tableau[row_idx, col_idx]
    end

    return ray
end
