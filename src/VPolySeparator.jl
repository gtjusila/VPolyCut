import SCIP

mutable struct VPolySeparator <: SCIP.AbstractSeparator
    called::Int64
    scipd::SCIP.SCIPData 
    row_num::Int64
    col_num::Int64
    lp_rows::Vector{Ptr{SCIP.SCIP_Row}}
    lp_cols::Vector{Ptr{SCIP.SCIP_Col}}
    VPolySeparator(inner) = new(0,inner,0,0,[],[])
end