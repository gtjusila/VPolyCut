import SCIP
mutable struct Ray
    coefficients::Vector{Ptr{SCIP.SCIP_Var}}
end