import SCIP

mutable struct VPolySeparator <: SCIP.AbstractSeparator
    called::Int64
    scipd::SCIP.SCIPData 
    VPolySeparator(inner) = new(0,inner)
end