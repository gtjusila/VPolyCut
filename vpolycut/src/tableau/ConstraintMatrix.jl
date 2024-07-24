# The ConstraintMatrix class allows us to store the constraint matrix of the LP at a particular 
# point of the solve. The Constraint Matrix stores the problem in SCIP Generatl Form, it
# distinguishes between the original variables and the slack variables.
using SCIP
using SparseArrays

mutable struct ConstraintMatrix
    data::SparseMatrixCSC{SCIP.SCIP_Real,Int}
    constants::SparseVector{SCIP.SCIP_Real,Int}
end

function ConstraintMatrix(nrows::Int, ncols::Int)
    return ConstraintMatrix(spzeros(nrows, ncols), spzeros(nrows))
end

function set_entry!(matrix::ConstraintMatrix, row::Int, col::Int, value::SCIP.SCIP_Real)
    matrix.data[row, col] = value
end

function set_constant!(matrix::ConstraintMatrix, row::Int, value::SCIP.SCIP_Real)
    matrix.constants[row] = value
end

function get_entry(matrix::ConstraintMatrix, row::Int, col::Int)::SCIP.SCIP_Real
    return matrix.data[row, col]
end

function get_constant(matrix::ConstraintMatrix, row::Int)::SCIP.SCIP_Real
    return matrix.constants[row]
end

function isacolumn(var::Variable)
    return isa(var, LPColumn)
end