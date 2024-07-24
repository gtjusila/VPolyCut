using SCIP
using VPolyCut
#
# Gomory Separator
#
struct GomorySeparator <: SCIP.AbstractSeparator
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

function include_separators(scip::SCIP.SCIPData, separator::Type{VPolyCut.IntersectionSeparator})
    VPolyCut.include_intersection_sepa(scip)
end