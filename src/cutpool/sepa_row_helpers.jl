using SCIP

# Add Row alpha * x <= beta to scip 
function add_sepa_row!(
    scip::SCIP.SCIPData,
    sepa::T,
    alpha::Vector{Float64},
    xs::Vector{Ptr{SCIP.SCIP_Var}},
    beta::Float64;
    valid_globaly::Bool = true,
    modifiable::Bool = false,
    removable::Bool = false
) where {T<:SCIP.AbstractSeparator}
    @assert length(alpha) == length(xs)

    new_row = Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRowSepa(
        scip,
        new_row,
        scip.sepas[sepa],
        "",
        -SCIP.SCIPinfinity(scip),
        beta,
        !valid_globaly,
        modifiable,
        removable
    )

    infeasible = Ref{SCIP.SCIP_Bool}(0)
    for (var, sol) in zip(xs, alpha)
        @assert !is_infinity(abs(sol))
        if is_non_zero(sol)
            SCIP.@SCIP_CALL SCIP.SCIPaddVarToRow(scip, new_row[], var, sol)
        end
    end

    #SCIP.@SCIP_CALL SCIP.SCIPprintRow(scip, new_row[], C_NULL)
    #SCIP.@SCIP_CALL SCIP.SCIPaddRow(scip, new_row[], true, infeasible)
    SCIP.@SCIP_CALL SCIP.SCIPaddPoolCut(scip, new_row[])
    #SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip, new_row)
    return new_row
end