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

    return Point(solution, 0.0)
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
            ray[basic_col] = value
        end
    end

    return ray
end

"""
    ConstructCornerPolyhedronContext

Context to construct a corner polyhedron
"""
struct ConstructCornerPolyhedronContext
    scip::SCIP.SCIPData
    basis_status::Vector{SCIP.SCIP_BASESTAT}
    nb_space::NonBasicSpace
    n_rows::Int64
    n_cols::Int64
    origin::Vector{Int64}
    complement::Vector{Int64}
    target::Vector{Int64}
end

"""
    CornerPolyhedron(scip::SCIP.SCIPData, nb_space::NonBasicSpace)::CornerPolyhedron

Create a corner polyhedron based on the LP currently loaded in SCIP projected to nb_space.
"""
function CornerPolyhedron(scip::SCIP.SCIPData, nb_space::NonBasicSpace)::CornerPolyhedron
    sol = get_solution_vector(scip)
    idx = 0 # initialize index to 0
    rays = Ray[]

    # Setup context
    basis_status = get_basis_status(scip)
    basic_indices = get_basis_indices(scip; sorted = false)
    n_cols = SCIP.SCIPgetNLPCols(scip)
    n_rows = SCIP.SCIPgetNLPRows(scip)
    origin = findall(x -> insorted(x, nb_space.nonbasic_indices), basic_indices)
    target = map(origin) do i
        id = basic_indices[i]
        return searchsortedfirst(nb_space.nonbasic_indices, id)
    end
    complement = map(origin) do i
        id = basic_indices[i]
        return id in nb_space.complemented_columns ? -1.0 : 1.0
    end
    context = ConstructCornerPolyhedronContext(
        scip, basis_status, nb_space, n_rows, n_cols, origin, complement, target
    )

    # Generate rays from nonbasic columns
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

"""
    construct_non_basic_ray_from_column(
    context::ConstructCornerPolyhedronContext,
    col::Ptr{SCIP.SCIP_Col},
    idx::Int64
)::Union{Ray,Nothing}

Construct a non basic ray from the column `col` at index `idx`
"""
function construct_non_basic_ray_from_column(
    context::ConstructCornerPolyhedronContext,
    col::Ptr{SCIP.SCIP_Col},
    idx::Int64
)::Union{Ray,Nothing}
    # It may be the case that column is nonbasic but we cannot move in any direction
    if is_EQ(SCIP.SCIPcolGetLb(col), SCIP.SCIPcolGetUb(col))
        return nothing
    end

    # the direction which we can move in
    direction = context.basis_status[idx] == SCIP.SCIP_BASESTAT_UPPER ? -1.0 : 1.0
    # Construct the ray
    nb_space_dim = length(context.nb_space.nonbasic_indices)
    ray = Ray(nb_space_dim, SCIP.SCIPgetColRedcost(context.scip, col))

    # First, the nonbasic column that is being considered should have entry `direction` in the ray.
    # We only care about is however if the column will be projected
    if insorted(idx, context.nb_space.nonbasic_indices)
        # get the index of idx under the projection
        i = searchsortedfirst(context.nb_space.nonbasic_indices, idx)
        # if the entry should be complemented then assign appropriate coefficient
        c = (idx in context.nb_space.complemented_columns) ? -1 : 1
        ray[i] = c * direction
    end

    # Get the actual data from the tableau column
    buffer_values = zeros(SCIP.SCIP_Real, context.n_rows)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvACol(
        context.scip, Cint(idx - 1), buffer_values, C_NULL, C_NULL)

    for (i, c, j) in zip(context.origin, context.complement, context.target)
        # we are looping through the rows i.e. through each basic variable
        # i is the indice of the row in the tableau, buffer[i] is the indice of the basic variable
        # c is either -1 or 1 depending on whether the column is complemented or not
        # j is the indice of the basic variable in the projection
        value = -direction * buffer_values[i]
        if !is_zero(value)
            ray[j] = c * value
        end
    end
    return ray
end

function construct_non_basic_ray_from_row(
    context::ConstructCornerPolyhedronContext,
    row::Ptr{SCIP.SCIP_Row},
    idx::Int64
)::Union{Ray,Nothing}
    if is_EQ(SCIP.SCIProwGetLhs(row), SCIP.SCIProwGetRhs(row))
        return nothing
    end

    # the direction which we can move in
    direction = context.basis_status[idx] == SCIP.SCIP_BASESTAT_UPPER ? -1.0 : 1.0
    # Construct the ray
    nb_space_dim = length(context.nb_space.nonbasic_indices)
    ray = Ray(nb_space_dim, -SCIP.SCIProwGetDualsol(row))

    if insorted(idx, context.nb_space.nonbasic_indices)
        # get the index of idx under the projection
        i = searchsortedfirst(context.nb_space.nonbasic_indices, idx)
        # if the entry should be complemented then assign appropriate coefficient
        c = (idx in context.nb_space.complemented_columns) ? -1 : 1
        ray[i] = c * direction
    end

    buffer_values = zeros(SCIP.SCIP_Real, context.n_rows)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvCol(
        context.scip, Cint(idx - context.n_cols - 1), buffer_values, C_NULL, C_NULL)

    for (i, c, j) in zip(context.origin, context.complement, context.target)
        # Same as in the column case
        value = -direction * buffer_values[i]
        if !is_zero(value)
            ray[j] = c * value
        end
    end

    return ray
end