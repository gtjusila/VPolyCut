# LPRow is a subtype of Variable
# It is used to represent slack variables arising from constraints
using SparseArrays

@kwdef mutable struct LPRow <: Variable
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    rhs::SCIP.SCIP_Real = -1
    lhs::SCIP.SCIP_Real = -1
    obj::SCIP.SCIP_Real = 0
    minus_row_activity::SCIP.SCIP_Real = -1
end

# upper and lower bounds for the slack variables are the left and right hand side of the constraint
get_ub(row::LPRow) = row.rhs
get_lb(row::LPRow) = row.lhs
get_sol(row::LPRow) = row.minus_row_activity # the solution of the slack variable is minus the row activity
get_obj(row::LPRow) = row.obj
get_basis_status(row::LPRow) = row.basis_status
get_symbolic_representation(row::LPRow) = :ROW
get_scip_index(row::LPRow) = row.scip_index
isarow(var::Variable) = isa(var, LPRow)

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

function set_obj!(row::LPRow, obj::SCIP.SCIP_Real)
    row.obj = obj
end

function set_sol!(row::LPRow, minus_row_activity::SCIP.SCIP_Real)
    row.minus_row_activity = minus_row_activity
end