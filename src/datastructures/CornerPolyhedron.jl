using SCIP
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
    complement::Vector{Int64}
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

    complement = map(relevant_indices) do i
        id = basic_indices[i]
        return ((id in nb_space.complemented_columns) ? -1 : 1)
    end

    n_rows::Int64 = SCIP.SCIPgetNLPRows(scip)
    n_cols::Int64 = SCIP.SCIPgetNLPCols(scip)
    dim = n_rows + n_cols
    context = ConstructCornerPolyhedronContext(
        scip, basis_status, n_rows, n_cols, dim, nb_space,
        length(nb_space.nonbasic_indices),
        complement, relevant_indices, to_put
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
    ray = Ray(context.nb_space_dim, SCIP.SCIPgetColRedcost(context.scip, col))
    if insorted(idx, context.nb_space.nonbasic_indices)
        i = searchsortedfirst(context.nb_space.nonbasic_indices, idx)
        ray[i] = ((idx in context.nb_space.complemented_columns) ? -1 : 1) * direction
    end

    buffer_values = zeros(SCIP.SCIP_Real, context.n_rows)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvACol(
        context.scip, Cint(idx - 1), buffer_values, C_NULL, C_NULL)
    for (i, c, j) in zip(context.relevant_indices, context.complement, context.to_put)
        value = -direction * buffer_values[i]
        if !is_zero(value)
            ray[j] = c * value
        end
    end
    return ray
end

function construct_non_basic_ray_from_row(
    context::ConstructCornerPolyhedronContext, row::Ptr{SCIP.SCIP_Row}, idx::Int64
)::Union{Ray,Nothing}
    if is_EQ(SCIP.SCIProwGetLhs(row), SCIP.SCIProwGetRhs(row))
        return nothing
    end

    direction = context.basis_status[idx] == SCIP.SCIP_BASESTAT_UPPER ? -1.0 : 1.0

    ray = Ray(context.nb_space_dim, -SCIP.SCIProwGetDualsol(row))
    if insorted(idx, context.nb_space.nonbasic_indices)
        i = searchsortedfirst(context.nb_space.nonbasic_indices, idx)
        ray[i] = ((idx in context.nb_space.complemented_columns) ? -1 : 1) * direction
    end

    buffer_values = zeros(SCIP.SCIP_Real, context.n_rows)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvCol(
        context.scip, Cint(idx - context.n_cols - 1), buffer_values, C_NULL, C_NULL)

    for (i, c, j) in zip(context.relevant_indices, context.complement, context.to_put)
        value = -direction * buffer_values[i]
        if !is_zero(value)
            ray[j] = c * value
        end
    end
    return ray
end