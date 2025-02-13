using SCIP
using DelimitedFiles
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

struct ConstructCornerPolyhedronContext
    scip::SCIP.SCIPData
    basis_status::Vector{SCIP.SCIP_BASESTAT}
    n_rows::Int64
    n_cols::Int64
    dim::Int64
    nb_space::NonBasicSpace
    nb_space_dim::Int64
    relevant_indices::Vector{Int64}
    to_put::Vector{Int64}
end
function construct_corner_polyhedron(
    scip::SCIP.SCIPData, nb_space::NonBasicSpace
)::CornerPolyhedron
    sol = get_solution_vector(scip)

    # Setup context for ray generation
    idx = 0
    basis_status = get_basis_status(scip)
    basic_indices = get_basis_indices(scip; sorted = false)

    # We do some efficiency here
    # Step 1 we need to find all elements of 
    relevant_indices = findall(x -> insorted(x, nb_space.nonbasic_indices), basic_indices)
    to_put = map(relevant_indices) do i
        id = basic_indices[i]
        return searchsortedfirst(nb_space.nonbasic_indices, id)
    end

    n_rows::Int64 = SCIP.SCIPgetNLPRows(scip)
    n_cols::Int64 = SCIP.SCIPgetNLPCols(scip)
    dim = n_rows + n_cols
    context = ConstructCornerPolyhedronContext(
        scip, basis_status, n_rows, n_cols, dim, nb_space,
        length(nb_space.nonbasic_indices),
        relevant_indices, to_put
    )

    rays = Vector{Ray}(undef, 0)

    cols = SCIP.SCIPgetLPCols(scip)
    cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, cols, n_cols)
    for col in cols
        idx += 1
        if basis_status[idx] == SCIP.SCIP_BASESTAT_BASIC
            continue
        end
        if basis_status[idx] == SCIP.SCIP_BASESTAT_ZERO
            throw(BasestatZeroEncountered())
        end
        ray = construct_non_basic_ray_from_column(context, col, idx)
        if !isnothing(ray)
            push!(rays, ray)
        end
    end

    rows = SCIP.SCIPgetLPRows(scip)
    rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, rows, n_rows)
    for row in rows
        idx += 1
        if basis_status[idx] == SCIP.SCIP_BASESTAT_BASIC
            continue
        end
        if basis_status[idx] == SCIP.SCIP_BASESTAT_ZERO
            throw(BasestatZeroEncountered())
        end
        ray = construct_non_basic_ray_from_row(context, row, idx)
        if !isnothing(ray)
            push!(rays, ray)
        end
    end
    return CornerPolyhedron(sol, rays)
end

function construct_non_basic_ray_from_column(
    context::ConstructCornerPolyhedronContext, col::Ptr{SCIP.SCIP_Col}, idx::Int64
)::Union{Ray,Nothing}
    if is_EQ(SCIP.SCIPcolGetLb(col), SCIP.SCIPcolGetUb(col))
        return nothing
    end

    direction = context.basis_status[idx] == SCIP.SCIP_BASESTAT_UPPER ? -1.0 : 1.0

    ray = Ray(context.nb_space_dim)
    if insorted(idx, context.nb_space.nonbasic_indices)
        i = searchsortedfirst(context.nb_space.nonbasic_indices, idx)
        ray[i] = direction
    end

    buffer_values = zeros(SCIP.SCIP_Real, context.n_rows)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvACol(
        context.scip, Cint(idx - 1), buffer_values, C_NULL, C_NULL)
    for (i, j) in zip(context.relevant_indices, context.to_put)
        value = -direction * buffer_values[i]
        if !is_zero(value)
            ray[j] = value
        end
    end
end

function construct_non_basic_ray_from_row(
    context::ConstructCornerPolyhedronContext, row::Ptr{SCIP.SCIP_Row}, idx::Int64
)::Union{Ray,Nothing}
    if is_EQ(SCIP.SCIProwGetLhs(row), SCIP.SCIProwGetRhs(row))
        return nothing
    end

    direction = context.basis_status[idx] == SCIP.SCIP_BASESTAT_UPPER ? -1.0 : 1.0

    ray = Ray(context.nb_space_dim)
    if insorted(idx, context.nb_space.nonbasic_indices)
        i = searchsortedfirst(context.nb_space.nonbasic_indices, idx)
        ray[i] = direction
    end

    buffer_values = zeros(SCIP.SCIP_Real, context.n_rows)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvCol(
        context.scip, Cint(idx - context.n_cols - 1), buffer_values, C_NULL, C_NULL)

    for (i, j) in zip(context.relevant_indices, context.to_put)
        value = -direction * buffer_values[i]
        if !is_zero(value)
            ray[j] = value
        end
    end
    return ray
end

"""
Create a Corner Polyhedron in the dimension of the problem in the SCIP Standard Form.
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
            println("Constructing Ray for: $(i)")
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
    ray = Ray(dim, var)

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
        if !is_EQ(tableau[row_idx, basic_col], 1.0)
            throw(AssumptionViolated())
        end
        value = -direction * tableau[row_idx, col_idx] * tableau[row_idx, basic_col]

        if !is_zero(value)
            println(basic_col)
            ray[basic_col] = value
        end
    end
    return ray
end
