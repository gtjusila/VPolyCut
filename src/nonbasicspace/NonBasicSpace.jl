
"""
    NonBasicSpace

A NonBasicSpace object represent the non-basic space associated with an LP tableau.
It stores information about the constraint matrix, the non-basic variables, the origin point of the non basic space and which columns are complemented
"""
struct NonBasicSpace
    nonbasic_indices::Vector{Int64}
    complemented_columns::Set{Int64}
    origin_point::Point
    constraint_matrix::ConstraintMatrix
    variable_pointers::Vector{Ptr{SCIP.SCIP_VAR}}
    auxiliary_objective::Vector{SCIP.SCIP_Real}
end

"""
    dimension(nbspace::NonBasicSpace)

Return the dimension of the non-basic space
"""
function dimension(nbspace::NonBasicSpace)
    return length(nbspace.nonbasic_indices)
end

"""
    NonBasicSpace(scip::SCIP.SCIPData)

Create a NonBasicSpace object from a given SCIP object
"""
function NonBasicSpace(scip::SCIP.SCIPData)::NonBasicSpace
    @info "Constructing NonBasicSpace"

    basis_status = get_basis_status(scip)
    nonbasic_indices = findall(
        x -> (x == SCIP.SCIP_BASESTAT_LOWER || x == SCIP.SCIP_BASESTAT_UPPER), basis_status
    )
    complemented_columns = findall(x -> x == SCIP.SCIP_BASESTAT_UPPER, basis_status)
    origin_point = get_solution_vector(scip)
    constraint_matrix = ConstraintMatrix(scip)
    variable_pointers = collect_variable_pointers(scip)
    auxiliary_objective = get_auxiliary_objective(scip)

    @info "Finished Constructing NonBasicSpace"
    return NonBasicSpace(
        nonbasic_indices,
        Set(complemented_columns),
        origin_point,
        constraint_matrix,
        variable_pointers,
        auxiliary_objective
    )
end

"""
    project_point_to_nonbasic_space(nbspace::NonBasicSpace, point::Point)

Project a point to the non-basic space
"""
function project_point_to_nonbasic_space(
    nbspace::NonBasicSpace,
    point::Point
)::Point
    # To project a point we need to remove the basic variables and complement the necessary non-basic columns
    point = substract(point, nbspace.origin_point)
    new_point = Point(dimension(nbspace), get_objective_value(point))
    for i in 1:dimension(nbspace)
        idx = nbspace.nonbasic_indices[i]
        new_point[i] = ((idx in nbspace.complemented_columns) ? -1 : 1) * point[idx]
    end
    return new_point
end

"""
    revert_cut_vector_to_nonbasic_space(nbspace::NonBasicSpace, cut_vector::Vector{SCIP.SCIP_Real})

Revert a cut vector to the original space (i.e. add the basic variables and undo the complementation) 
"""
function revert_cut_vector_to_original_space(
    nbspace::NonBasicSpace,
    cut_vector::Vector{SCIP.SCIP_Real}
)::Vector{SCIP.SCIP_Real}
    new_cut_vector = zeros(SCIP.SCIP_Real, length(nbspace.origin_point))

    # Add the basic variables
    for (idx, val) in enumerate(cut_vector)
        new_cut_vector[nbspace.nonbasic_indices[idx]] = val
    end

    # Undo the complementation
    for idx in nbspace.complemented_columns
        new_cut_vector[idx] = -new_cut_vector[idx]
    end

    return new_cut_vector
end

"""
    get_auxiliary_objective(scip::SCIP.SCIPData)

Give an objective such that the current SCIP basic feasible solution is the unique minimizer over the polyhedron for this objective
"""
function get_auxiliary_objective(scip::SCIP.SCIPData)
    # Formula from deepseek :)
    # Should eventually be proven but seems correct
    nlpcols = Int64(SCIP.SCIPgetNLPCols(scip))
    lp_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, lp_cols, nlpcols)

    objective = zeros(nlpcols)
    var_ptrs = []
    statuses = []
    for (idx, col) in enumerate(lp_cols)
        push!(var_ptrs, SCIP.SCIPcolGetVar(col))
        push!(statuses, SCIP.SCIPcolGetBasisStatus(col))
        if SCIP.SCIPcolGetBasisStatus(col) == SCIP.SCIP_BASESTAT_LOWER
            objective[idx] = 1
        elseif SCIP.SCIPcolGetBasisStatus(col) == SCIP.SCIP_BASESTAT_UPPER
            objective[idx] = -1
        end

        nnonz = SCIP.SCIPcolGetNNonz(col)
        nonzero_rows = SCIP.SCIPcolGetRows(col)
        nonzero_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, nonzero_rows, nnonz)
        nonzero_vals = SCIP.SCIPcolGetVals(col)
        nonzero_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, nonzero_vals, nnonz)
        for (row_idx, row) in enumerate(nonzero_rows)
            if SCIP.SCIProwGetBasisStatus(row) == SCIP.SCIP_BASESTAT_UPPER
                objective[idx] += -nonzero_vals[row_idx]
            elseif SCIP.SCIProwGetBasisStatus(row) == SCIP.SCIP_BASESTAT_LOWER
                objective[idx] += nonzero_vals[row_idx]
            end
        end
    end

    # Verify this
    solution1 = get_solution_vector(scip)
    objvalue = dot(solution1.coordinates[1:length(objective)], objective)

    @info "Objective $(objvalue)"

    SCIP.@SCIP_CALL SCIP.SCIPstartDive(scip)
    for (idx, var) in enumerate(var_ptrs)
        SCIP.@SCIP_CALL SCIP.SCIPchgVarObjDive(scip, var, objective[idx])
    end
    lperror = Ref{SCIP.SCIP_Bool}(0)
    cutoff = Ref{SCIP.SCIP_Bool}(0)
    SCIP.SCIPsolveDiveLP(scip, -1, lperror, cutoff)
    @info "Obj $(SCIP.SCIPgetLPObjval(scip))"
    solution2 = get_solution_vector(scip)
    nlpcols = Int64(SCIP.SCIPgetNLPCols(scip))
    lp_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, lp_cols, nlpcols)
    statuses2 = []
    for col in lp_cols
        push!(statuses2, SCIP.SCIPcolGetBasisStatus(col))
    end

    for i in 1:length(statuses)
        if (statuses[i] != statuses2[i])
            @error "i = $i. Nonmatching status $(statuses[i]) not equal $(statuses2[i])"
            throw(LPError())
        else
            @info "Status match. $(solution1[i]) $(solution2[i])"
        end
    end
    #@assert solution1.coordinates == solution2.coordinates
    @assert equal(solution1.coordinates, solution2.coordinates; comp = is_EQ)
    SCIP.@SCIP_CALL SCIP.SCIPendDive(scip)

    return objective
end

function equal(v1::AbstractVector, v2::AbstractVector; comp::Function = (x, y) -> x == y)
    if length(v1) != length(v2)
        return false
    end
    for i in 1:length(v1)
        if !comp(v1[i], v2[i])
            @info "Differ on $(i)-th component"
            @info "$(v1[i]) vs $(v2[i])"
            return false
        end
    end
    return true
end
