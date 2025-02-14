using SCIP

### Basis Related Utilities ###
"""
    get_basis_status(scip::SCIP.SCIPData)

Return vector of basis status for all the variables and rows in the current LP solution.
We reverse the basis status for the rows to follow the convention that column with basis status upper should be complemented
"""
function get_basis_status(scip::SCIP.SCIPData)
    n_rows::Int64 = SCIP.SCIPgetNLPRows(scip)
    n_cols::Int64 = SCIP.SCIPgetNLPCols(scip)
    dim = n_rows + n_cols
    basis_status = Vector{SCIP.SCIP_BASESTAT}(undef, dim)
    idx = 1

    cols = SCIP.SCIPgetLPCols(scip)
    cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, cols, n_cols)
    for col in cols
        basis_status[idx] = SCIP.SCIPcolGetBasisStatus(col)
        idx += 1
    end

    rows = SCIP.SCIPgetLPRows(scip)
    rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, rows, n_rows)
    for row in rows
        status = SCIP.SCIProwGetBasisStatus(row)
        if status == SCIP.SCIP_BASESTAT_LOWER
            status = SCIP.SCIP_BASESTAT_UPPER
        elseif status == SCIP.SCIP_BASESTAT_UPPER
            status = SCIP.SCIP_BASESTAT_LOWER
        end
        basis_status[idx] = status
        idx += 1
    end

    return basis_status
end

"""
    get_basis_indices(scip::SCIP.SCIPData; sorted = true)

Return the basis indices for the current LP solution. The indices are 1-based and sorted by default.
"""
function get_basis_indices(scip::SCIP.SCIPData; sorted = true)::Vector{Int64}
    n_rows::Int64 = SCIP.SCIPgetNLPRows(scip)
    n_cols::Int64 = SCIP.SCIPgetNLPCols(scip)
    nonbasic_indices = zeros(Cint, n_rows)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, pointer(nonbasic_indices))
    nonbasic_indices = map(nonbasic_indices) do x
        if x >= 0
            x += 1
        else
            x *= -1
            x += n_cols
        end
        return x
    end
    # Finally convert indexing to Int64
    if sorted
        sort!(nonbasic_indices)
    end

    return nonbasic_indices
end

### LP Solution Methods ###

"""
    get_solution_vector(scip::SCIP.SCIPData)::Point

Given a SCIPData object, return a Point object representing the solution vector.
Our solution is (x, -y) where x is the solution to the problem variables and y is the row activity.
"""
function get_solution_vector(scip::SCIP.SCIPData)::Point
    n_rows::Int64 = SCIP.SCIPgetNLPRows(scip)
    n_cols::Int64 = SCIP.SCIPgetNLPCols(scip)
    dim = n_rows + n_cols
    sol = Point(dim, SCIP.SCIPgetSolOrigObj(scip, C_NULL))
    idx = 1

    cols = SCIP.SCIPgetLPCols(scip)
    cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, cols, n_cols)
    for col in cols
        val = SCIP.SCIPcolGetPrimsol(col)
        if !is_zero(val)
            sol[idx] = val
        end
        idx += 1
    end

    rows = SCIP.SCIPgetLPRows(scip)
    rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, rows, n_rows)
    for row in rows
        val = -SCIP.SCIPgetRowActivity(scip, row)
        if !is_zero(val)
            sol[idx] = val
        end
        idx += 1
    end

    return sol
end
