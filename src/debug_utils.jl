include("VPolySeparator.jl")
# [CLEAN] Debug Function add_probing_bound_via_row
function add_probing_bound_via_row(sepa::VPolySeparator, var::Ptr{SCIP.SCIP_Var}; upper = nothing, lower = nothing)
    if !isnothing(upper)
        row = Ref{Ptr{SCIP.SCIP_Row}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRow(sepa.scipd, row, "", -SCIP.SCIPinfinity(sepa.scipd),upper, true, false, true)
        SCIP.@SCIP_CALL SCIP.SCIPaddVarToRow(sepa.scipd, row[],var,1)
        SCIP.@SCIP_CALL SCIP.SCIPaddRowProbing(sepa.scipd, row[])
    end
    if !isnothing(lower)
        row = Ref{Ptr{SCIP.SCIP_Row}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRow(sepa.scipd, row, "", lower ,SCIP.SCIPinfinity(sepa.scipd), true, false, true)
        SCIP.@SCIP_CALL SCIP.SCIPaddVarToRow(sepa.scipd, row[],var,1)
        SCIP.@SCIP_CALL SCIP.SCIPaddRowProbing(sepa.scipd, row[])
    end
end