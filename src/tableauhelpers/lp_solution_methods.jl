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

    sol = Point(dim)

    cols = SCIP.SCIPgetLPCols(scip)
    cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, cols, n_cols)
    for i in 1:n_cols
        col = cols[i]
        sol[i] = SCIP.SCIPcolGetPrimsol(col)
    end

    rows = SCIP.SCIPgetLPRows(scip)
    rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, rows, n_rows)
    for i in 1:n_rows
        row = rows[i]
        sol[n_cols + i] = -SCIP.SCIPgetRowActivity(scip, row)
    end

    return sol
end