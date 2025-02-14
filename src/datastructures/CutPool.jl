using SCIP
using Lazy
"""
A Cut of the form `coefficients * x <= rhs`
"""
@kwdef mutable struct Cut
    coefficients::Vector{SCIP.SCIP_Real} = []
    rhs::SCIP.SCIP_Real = []
end

function get_coefficients(cut::Cut)
    return cut.coefficients
end

function get_rhs(cut::Cut)
    return cut.rhs
end

@kwdef mutable struct CutPool <: AbstractVector{Cut}
    problem_var_pointers::Vector{Ptr{SCIP.SCIP_Var}} = []
    cuts::Vector{Cut} = []
end

@forward CutPool.cuts Base.push!, Base.size

function add_all_cuts!(
    scip::SCIP.SCIPData, cutpool::CutPool, sepa::T
) where {T<:SCIP.AbstractSeparator}
    for cut in cutpool.cuts
        row = add_sepa_row!(
            scip,
            sepa,
            get_coefficients(cut),
            cutpool.problem_var_pointers,
            get_rhs(cut)
        )
        SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip, row)
    end
end