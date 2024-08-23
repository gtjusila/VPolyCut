using SCIP
using VPolyhedralCut

#
# Gomory Separator
#
struct GomorySeparator <: SCIP.AbstractSeparator
end

function include_separators(
    scip::SCIP.SCIPData,
    separator::Type{T}
) where {T<:SCIP.AbstractSeparator}
    error("include_separators not supported for type $T")
end

function include_separators(scip::SCIP.SCIPData, separator::Type{GomorySeparator})
    setter = (par, val) -> SCIP.set_parameter(scip, par, val)
    setter("separating/gmi/freq", 0)
    setter("separating/gmi/priority", 9999)
    setter("separating/gmi/maxsuppabs", 5000)
    setter("separating/gmi/dynamiccuts", false)
    setter("separating/gmi/maxsupprel", 1.0)
    setter("separating/gmi/forcecuts", true)
end

function include_separators(
    scip::SCIP.SCIPData, separator::Type{VPolyhedralCut.IntersectionSeparator}
)
    VPolyhedralCut.include_intersection_sepa(scip)
end

function include_separators(
    scip::SCIP.SCIPData, separator::Type{VPolyhedralCut.VPolyhedralSeparator}
)
    VPolyhedralCut.include_vpolyhedral_sepa(scip; n_leaves=64)
end