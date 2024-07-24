# LPRow is a subtype of Variable
# It is used to represent slack variables arising from constraints

@kwdef mutable struct LPRow <: Variable
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    rhs::SCIP.SCIP_Real = -1
    lhs::SCIP.SCIP_Real = -1
    obj_pointer::Ptr{SCIP.SCIP_ROW} = C_NULL
    coefficient::Vector{SCIP.SCIP_Real} = []
    constant::SCIP.SCIP_Real = 0
    row_activity::SCIP.SCIP_Real = -1
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
    return -row.row_activity
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

function get_row_coefficients(row::LPRow)::Vector{SCIP.SCIP_Real}
    return row.coefficient
end

function get_constant(row::LPRow)::SCIP.SCIP_Real
    return row.constant
end
"""
Populate Row with information from SCIP, Tableau is necessary to get column information
"""
function populate!(scip::SCIP.SCIPData, scip_ptr::Ptr{SCIP.SCIP_ROW}, row::LPRow)
    row.scip_index = SCIP.SCIProwGetIndex(scip_ptr)
    row.basis_status = SCIP.SCIProwGetBasisStatus(scip_ptr)
    row.rhs = SCIP.SCIProwGetRhs(scip_ptr)
    row.lhs = SCIP.SCIProwGetLhs(scip_ptr)
    row.row_activity = SCIP.SCIPgetRowLPActivity(scip, scip_ptr)
    row.constant = SCIP.SCIProwGetConstant(scip_ptr)
end

function set_row_coefficient!(row::LPRow, coefficient::Vector{SCIP.SCIP_Real})
    row.coefficient = coefficient
end

