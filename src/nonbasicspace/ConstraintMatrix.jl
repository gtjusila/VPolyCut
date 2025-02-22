# The ConstraintMatrix class allows us to store the constraint matrix of the LP at a particular 
# point of the solve. The Constraint Matrix stores the problem in SCIP Generatl Form, it
# distinguishes between the original variables and the slack variables.
using SCIP
using SparseArrays

struct ConstraintMatrix <: AbstractMatrix{SCIP.SCIP_Real}
    data::SparseMatrixCSC{SCIP.SCIP_Real,Int}
    constants::SparseVector{SCIP.SCIP_Real,Int}
end

get_entry(matrix::ConstraintMatrix, row::Int, col::Int) = matrix.data[row, col]
get_constant(matrix::ConstraintMatrix, row::Int) = matrix.constants[row]

function ConstraintMatrix(nrows::Int, ncols::Int)
    return ConstraintMatrix(spzeros(nrows, ncols), spzeros(nrows))
end

"""
    ConstraintMatrix(scip::SCIP.SCIPData)::ConstraintMatrix

Create a constraint matrix based on the LP currently loaded in SCIP.
"""
function ConstraintMatrix(scip::SCIP.SCIPData)::ConstraintMatrix
    n_rows = SCIP.SCIPgetNLPRows(scip)
    n_cols = SCIP.SCIPgetNLPCols(scip)
    constants = zeros(n_rows)

    # We loop row by row
    rows = SCIP.SCIPgetLPRows(scip)
    rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, rows, n_rows)
    total = 0

    for (idx, row) in enumerate(rows)
        nnonz = SCIP.SCIProwGetNNonz(row)
        total += nnonz
    end

    sparse_row = Int64[]
    sparse_col = Int64[]
    sparse_val = SCIP.SCIP_Real[]
    sizehint!(sparse_row, total)
    sizehint!(sparse_col, total)
    sizehint!(sparse_val, total)

    for (idx, row) in enumerate(rows)
        nnonz = SCIP.SCIProwGetNNonz(row)

        nonzero_columns = SCIP.SCIProwGetCols(row)
        nonzero_columns = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, nonzero_columns, nnonz)
        nonzero_columns = map(x -> Int64(SCIP.SCIPcolGetLPPos(x)) + 1, nonzero_columns)

        values = SCIP.SCIProwGetVals(row)
        values = unsafe_wrap(Vector{SCIP.SCIP_Real}, values, nnonz)
        append!(sparse_row, ones(Int64, nnonz) * idx)
        append!(sparse_col, nonzero_columns)
        append!(sparse_val, values)
        # Set the constant term for the row 
        cons = SCIP.SCIProwGetConstant(row)
        constants[idx] = cons
    end
    data = sparse(sparse_row, sparse_col, sparse_val, n_rows, n_cols)
    return ConstraintMatrix(data, constants)
end

function set_entry!(matrix::ConstraintMatrix, row::Int, col::Int, value::SCIP.SCIP_Real)
    matrix.data[row, col] = value
end

function set_constant!(matrix::ConstraintMatrix, row::Int, value::SCIP.SCIP_Real)
    matrix.constants[row] = value
end
function Base.size(matrix::ConstraintMatrix)
    return size(matrix.data)
end