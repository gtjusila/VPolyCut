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

const global buffer = Ref{Vector{SCIP.SCIP_Real}}(Vector{SCIP.SCIP_Real}());
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

    # The origin vector contains all of the indices in our original space which will contain
    # relevant information for us to generate rays
    # The key insight here is these are exactly the nonbasic_indices in the original tbaleau that
    # have become basic in the current disjunction
    origin = findall(x -> insorted(x, nb_space.nonbasic_indices), basic_indices)
    # The target is which indices each of the indices in origin will be mapped to in the nonbasic space 
    target = map(origin) do i
        id = basic_indices[i]
        return searchsortedfirst(nb_space.nonbasic_indices, id)
    end
    # The complement vector contains information on whether or not the relevant indices have to be complemented or not
    complement = map(origin) do i
        id = basic_indices[i]
        return id in nb_space.complemented_columns ? -1.0 : 1.0
    end
    context = ConstructCornerPolyhedronContext(
        scip, basis_status, nb_space, n_rows, n_cols, origin, complement, target
    )
    if length(buffer[]) != context.n_rows
        # Global buffer is not set so set it
        buffer[] = zeros(SCIP.SCIP_Real, context.n_rows)
    end

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

    # Get the actual data from the tableau column
    global buffer
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvACol(
        context.scip, Cint(idx - 1), buffer[], C_NULL, C_NULL)

    # To make things efficient we must first determine which entry of the rays will be nonzero 
    non_zero = findall(i -> !is_zero(buffer[][i]), context.origin)
    indices = context.origin[non_zero]
    comp = context.complement[non_zero]
    target = context.target[non_zero]

    # Now we construct the ray
    nb_space_dim = length(context.nb_space.nonbasic_indices)
    ray = Ray(
        sparsevec(
            target,
            -direction * (comp .* buffer[][indices]),
            nb_space_dim
        ),
        SCIP.SCIPgetColRedcost(context.scip, col)
    )

    # Finally, the nonbasic column that is being considered should have entry `direction` in the ray.
    # We only care about is however if the column will be projected
    if insorted(idx, context.nb_space.nonbasic_indices)
        # get the index of idx under the projection
        i = searchsortedfirst(context.nb_space.nonbasic_indices, idx)
        # if the entry should be complemented then assign appropriate coefficient
        c = (idx in context.nb_space.complemented_columns) ? -1 : 1
        ray[i] = c * direction
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

    # Get the actual data from the tableau column
    global buffer
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvCol(
        context.scip, Cint(idx - context.n_cols - 1), buffer[], C_NULL, C_NULL)

    # Same as above
    non_zero = findall(i -> !is_zero(buffer[][i]), context.origin)
    indices = context.origin[non_zero]
    comp = context.complement[non_zero]
    target = context.target[non_zero]

    # Construct the ray
    nb_space_dim = length(context.nb_space.nonbasic_indices)
    ray = Ray(
        sparsevec(
            target,
            -direction * (comp .* buffer[][indices]),
            nb_space_dim
        ),
        -SCIP.SCIProwGetDualsol(row)
    )

    if insorted(idx, context.nb_space.nonbasic_indices)
        # get the index of idx under the projection
        i = searchsortedfirst(context.nb_space.nonbasic_indices, idx)
        # if the entry should be complemented then assign appropriate coefficient
        c = (idx in context.nb_space.complemented_columns) ? -1 : 1
        ray[i] = c * direction
    end

    return ray
end