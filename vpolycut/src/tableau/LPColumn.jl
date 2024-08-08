# LPColumn is a subtype of Variable
# It is used to represent the original problem variables

@kwdef mutable struct LPColumn <: Variable
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    ub::SCIP.SCIP_Real = -1
    lb::SCIP.SCIP_Real = -1
    sol::SCIP.SCIP_Real = -1
    obj::SCIP.SCIP_Real = -1
    variable_pointer::Ptr{SCIP.SCIP_VAR} = C_NULL
end

get_lb(col::LPColumn) = col.lb
get_ub(col::LPColumn) = col.ub
get_obj(col::LPColumn) = col.obj
get_sol(col::LPColumn) = col.sol   # The solution of a problem variable is the primal solution
get_basis_status(col::LPColumn) = col.basis_status
get_symbolic_representation(col::LPColumn) = :COLUMN
get_scip_index(col::LPColumn) = col.scip_index
get_var_pointer(col::LPColumn) = col.variable_pointer # Column Var Pointer would point to the SCIP variable associated with the column
isacolumn(var::Variable) = isa(var, LPColumn)
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

function set_obj!(col::LPColumn, obj::SCIP.SCIP_Real)
    col.obj = obj
end

function set_var_pointer!(col::LPColumn, pointer::Ptr{SCIP.SCIP_VAR})
    col.variable_pointer = pointer
end