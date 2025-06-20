using SCIP
function generate_split_cuts(
    sepa::VPCSeparator
)
    scipd = sepa.shared_data.scipd
    ncols = SCIP.SCIPgetNLPcols(scipd)
    nrows = SCIP.SCIPgetNLProws(scipd)
    cols_data = SCIP.SCIPgetLPCols(scipd)
    cols_data = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, cols_data, ncols)
    rows_data = SCIP.SCIPgetLPRows(scipd)
    rows_data = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, rows_data, nrows)
end