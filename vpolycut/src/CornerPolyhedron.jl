using SCIP

struct CornerPolyhedron
    lp_sol::Point
    lp_rays::Vector{Ray}
end

"""
Create a Corner Polyhedron in the dimension of the problem in the SCIP Standard Form
"""
function construct_corner_polyhedron(tableau::Tableau)::CornerPolyhedron
    # Initiate a vector to collect corner polyhedron ray
    rays = get_non_basic_rays(tableau)
    sol = get_solution_vector(tableau)

    return CornerPolyhedron(sol, rays)
end

function get_solution_vector(tableau::Tableau)::Point
    dim = get_nvars(tableau)
    solution = zeros(dim)

    for i in 1:dim
        var = get_var_from_column(tableau, i)
        solution[i] = get_sol(var)
    end

    return solution
end

function get_non_basic_rays(tableau::Tableau)::Vector{Ray}
    ray_collection = Vector{Ray}(undef, 0)

    for i in 1:get_nvars(tableau)
        var = get_var_from_column(tableau, i)
        if !is_basic(var)
            ray = construct_non_basic_ray(tableau, var)
            push!(ray_collection, ray)
        end
    end

    return ray_collection
end

"""
Construct non basic ray from the ith column
"""
function construct_non_basic_ray(tableau::Tableau, var::Variable)::Ray
    direction = 1.0

    if is_at_upper_bound(var)
        direction = -1.0
    elseif is_at_lower_bound(var)
        direction = 1.0
    else
        #Safekeeping: Should Never Happen unless polyhedron is not pointed
        error("Invalid basis status encountered: $(get_basis_status(var))")
    end

    # SCIP assumes that the constraint matrix is in the form [A I] where I
    # are columns corresponding to the slack variables. Hence, if the 
    # nonbasic column is a slack variable, then the direction is reversed
    col_idx = get_column_from_var(tableau, var)
    if col_idx > get_noriginalcols(tableau)
        direction = -direction
    end

    # Construct ray r
    dim = get_nvars(tableau)
    ray = zeros(dim)

    ray[col_idx] = direction

    #
    # Suppose the tableau is 
    # 1x1 + 1x2 + 1s = 0
    # where x1 is basic, x2 is at its lower bound and s is at
    # its upper bound. Then the ray corresponding to 
    # x2 is [-1 1 0] and the ray corresponding to s is [1 0 -1]
    for row_idx in 1:get_nbasis(tableau)
        basic_var = get_var_from_row(tableau, row_idx)
        basic_col = get_column_from_var(tableau, basic_var)
        value = -direction * tableau[row_idx, col_idx]
        ray[basic_col] = value
    end

    return ray
end
