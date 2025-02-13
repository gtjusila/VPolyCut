using SCIP

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