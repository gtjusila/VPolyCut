#This file contains the interface a Variable type
using SCIP

abstract type Variable end

function is_basic(var::Variable)::Bool
    return get_basis_status(var) == SCIP.SCIP_BASESTAT_BASIC
end

function is_at_upper_bound(var::Variable)::Bool
    return get_basis_status(var) == SCIP.SCIP_BASESTAT_UPPER
end

function is_at_lower_bound(var::Variable)::Bool
    return get_basis_status(var) == SCIP.SCIP_BASESTAT_LOWER
end

function is_zero(var::Variable)::Bool
    return get_basis_status(var) == SCIP.SCIP_BASESTAT_ZERO
end

# To implement are the following functions
function get_basis_status(var::Variable)::SCIP.SCIP_BASESTAT
    error("get_basis_status not implemented for type ", typeof(var))
end

function get_ub(var::Variable)::SCIP.SCIP_Real
    error("get_ub not implemented for type ", typeof(var))
end

function get_lb(var::Variable)::SCIP.SCIP_Real
    error("get_lb not implemented for type ", typeof(var))
end

function get_sol(var::Variable)::SCIP.SCIP_Real
    error("get_sol not implemented for type ", typeof(var))
end

"""
This function should return a unique representation for the subtype of Variable
"""
function get_symbolic_representation(var::Variable)::Symbol
    error("get_symbolic_representation not implemented for type ", typeof(var))
end

function get_scip_index(var::Variable)::Int
    error("get_scip_index not implemented for type ", typeof(var))
end

function populate!(scip::SCIP.SCIPData, scip_ptr::Ptr, var::Variable)
    error("populate! not implemented for type ", typeof(var))
end