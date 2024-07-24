# LPRow is a subtype of Variable
# It is used to represent slack variables arising from constraints
using SparseArrays

@kwdef mutable struct LPRow <: Variable
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    rhs::SCIP.SCIP_Real = -1
    lhs::SCIP.SCIP_Real = -1
    minus_row_activity::SCIP.SCIP_Real = -1
end

# upper and lower bounds for the slack variables are the left and right hand side of the constraint
function get_ub(row::LPRow)::SCIP.SCIP_Real
    return row.rhs
end

function get_lb(row::LPRow)::SCIP.SCIP_Real
    return row.lhs
end

# the solution of the slack variable is minus the row activity 
function get_sol(row::LPRow)::SCIP.SCIP_Real
    return row.minus_row_activity
end

function get_basis_status(row::LPRow)::SCIP.SCIP_BASESTAT
    return row.basis_status
end

function get_symbolic_representation(row::LPRow)::Symbol
    return :ROW
end

function get_scip_index(row::LPRow)::Int
    return row.scip_index
end

function set_scip_index!(row::LPRow, index::Integer)
    row.scip_index = index
end

function set_basis_status!(row::LPRow, status::SCIP.SCIP_BASESTAT)
    row.basis_status = status
end

function set_ub!(row::LPRow, ub::SCIP.SCIP_Real)
    row.rhs = ub
end

function set_lb!(row::LPRow, lb::SCIP.SCIP_Real)
    row.lhs = lb
end

function set_sol!(row::LPRow, minus_row_activity::SCIP.SCIP_Real)
    row.minus_row_activity = minus_row_activity
end

function isarow(var::Variable)
    return isa(var, LPRow)
end