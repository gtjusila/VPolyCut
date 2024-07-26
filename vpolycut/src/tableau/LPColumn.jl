# LPColumn is a subtype of Variable
# It is used to represent the original problem variables

@kwdef mutable struct LPColumn <: Variable
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    ub::SCIP.SCIP_Real = -1
    lb::SCIP.SCIP_Real = -1
    sol::SCIP.SCIP_Real = -1
    variable_pointer::Ptr{SCIP.SCIP_VAR} = C_NULL
end

function get_lb(col::LPColumn)::SCIP.SCIP_Real
    return col.lb
end

function get_ub(col::LPColumn)::SCIP.SCIP_Real
    return col.ub
end

# The solution of a problem variable is the primal solution
function get_sol(col::LPColumn)::SCIP.SCIP_Real
    return col.sol
end

function get_basis_status(col::LPColumn)::SCIP.SCIP_BASESTAT
    return col.basis_status
end

function get_symbolic_representation(row::LPColumn)::Symbol
    return :COLUMN
end

function get_scip_index(col::LPColumn)::Int
    return col.scip_index
end

"""
Column Var Pointer would point to the SCIP variable associated with the column
"""
function get_var_pointer(col::LPColumn)::Ptr{SCIP.SCIP_VAR}
    return col.variable_pointer
end

function set_scip_index!(col::LPColumn, index::Integer)
    col.scip_index = index
end

function set_basis_status!(col::LPColumn, status::SCIP.SCIP_BASESTAT)
    col.basis_status = status
end

function set_ub!(col::LPColumn, ub::SCIP.SCIP_Real)
    col.ub = ub
end

function set_lb!(col::LPColumn, lb::SCIP.SCIP_Real)
    col.lb = lb
end

function set_sol!(col::LPColumn, sol::SCIP.SCIP_Real)
    col.sol = sol
end

function set_var_pointer!(col::LPColumn, pointer::Ptr{SCIP.SCIP_VAR})
    col.variable_pointer = pointer
end
