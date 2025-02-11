using SCIP: SCIP

struct CornerPolyhedron
    lp_sol::Point
    lp_rays::Vector{Ray}
end

function get_lp_sol(corner::CornerPolyhedron)::Point
    return corner.lp_sol
end

function get_lp_rays(corner::CornerPolyhedron)::Vector{Ray}
    return corner.lp_rays
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

    return Point(solution)
end

function get_non_basic_rays(tableau::Tableau)::Vector{Ray}
    ray_collection = Vector{Ray}(undef, 0)

    for i in 1:get_nvars(tableau)
        var = get_var_from_column(tableau, i)
        if !is_basic(var)
            ray = construct_non_basic_ray(tableau, var)
            if !isnothing(ray)
                push!(ray_collection, ray)
            end
        end
    end
    return ray_collection
end

"""
Construct non basic ray from the ith column
"""
function construct_non_basic_ray(
    tableau::Tableau, var::Variable
)::Union{Ray,Nothing}
    direction = 1.0

    if is_EQ(get_ub(var), get_lb(var))
        return nothing
    end
    if is_at_upper_bound(var)
        direction = -1.0
    elseif is_at_lower_bound(var)
        direction = 1.0
    else
        #Safekeeping: Should Never Happen unless polyhedron is not pointed
        error("Invalid basis status encountered: $(get_basis_status(var))")
    end

    col_idx = get_column_from_var(tableau, var)
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
        # When we do point ray collection, we compliment the columns in according to the original tableau of the node
        # it may be the case that a basic column is complemented in this case the assumption
        # that the basic column coefficients is +1 is not valid, i.e. the last term in the following line may be -1
        value = -direction * tableau[row_idx, col_idx] * tableau[row_idx, basic_col]
        ray[basic_col] = value
    end

    return Ray(ray, var)
end
