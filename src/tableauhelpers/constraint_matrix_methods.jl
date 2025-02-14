using SCIP

function get_constraint_matrix(scip::SCIP.SCIPData)::ConstraintMatrix
    n_rows = SCIP.SCIPgetNLPRows(scip)
    n_cols = SCIP.SCIPgetNLPCols(scip)
    matrix = ConstraintMatrix(Int64(n_rows), Int64(n_cols))

    rows = SCIP.SCIPgetLPRows(scip)
    rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, rows, n_rows)

    for (row_index, row) in enumerate(rows)
        row_entries = fetch_non_zero_entries_from_row(row)
        for (col_index, entry) in row_entries
            set_entry!(matrix, row_index, col_index, entry)
        end
        set_constant!(matrix, row_index, SCIP.SCIProwGetConstant(row))
    end
    return matrix
end

function fetch_non_zero_entries_from_row(
    row_pointer::Ptr{SCIP.SCIP_ROW}
)::Vector{Tuple{Int,SCIP.SCIP_Real}}
    buffer::Vector{Tuple{Int,SCIP.SCIP_Real}} = []
    nnonzeros = SCIP.SCIProwGetNNonz(row_pointer)

    nonzero_columns = SCIP.SCIProwGetCols(row_pointer)
    nonzero_columns = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, nonzero_columns, nnonzeros)
    nonzero_entries = SCIP.SCIProwGetVals(row_pointer)
    nonzero_entries = unsafe_wrap(Vector{SCIP.SCIP_Real}, nonzero_entries, nnonzeros)

    for i in 1:nnonzeros
        col_index = SCIP.SCIPcolGetLPPos(nonzero_columns[i]) + 1
        col_entry = nonzero_entries[i]
        if col_index > 0
            push!(buffer, (col_index, col_entry))
        end
    end

    return buffer
end

function get_problem_variables_pointers(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_VAR}}
    n_cols = SCIP.SCIPgetNLPCols(scip)
    cols = SCIP.SCIPgetLPCols(scip)
    cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, cols, n_cols)
    pointers = Ptr{SCIP.SCIP_VAR}[]
    for col in cols
        ptr = SCIP.SCIPcolGetVar(col)
        push!(pointers, ptr)
    end
    return pointers
end