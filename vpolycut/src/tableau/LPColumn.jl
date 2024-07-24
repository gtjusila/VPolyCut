# LPColumn is a subtype of Variable
# It is used to represent the original problem variables

@kwdef mutable struct LPColumn <: Variable
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    ub::SCIP.SCIP_Real = -1
    lb::SCIP.SCIP_Real = -1
    sol::SCIP.SCIP_Real = -1
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
Populate column with information from SCIP
"""
function populate!(scip::SCIP.SCIPData, scip_ptr::Ptr{SCIP.SCIP_COL}, col::LPColumn)
    col.scip_index = SCIP.SCIPcolGetIndex(scip_ptr)
    col.basis_status = SCIP.SCIPcolGetBasisStatus(scip_ptr)
    col.ub = SCIP.SCIPcolGetUb(scip_ptr)
    col.lb = SCIP.SCIPcolGetLb(scip_ptr)
    col.sol = SCIP.SCIPcolGetPrimsol(scip_ptr)
end