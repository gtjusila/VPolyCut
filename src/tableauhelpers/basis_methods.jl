using SCIP

function get_non_basic_indices(scip::SCIP.SCIPData)::Vector{Int64}
    n_cols::Int64 = SCIP.SCIPgetNLPCols(scip)
    cols = SCIP.SCIPgetLPCols(scip)
    cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, cols, n_cols)
    non_basic_indices = Vector{Int64}()
    for i in 1:n_cols
        col = cols[i]
        if SCIP.SCIPcolIsBasic(col) == 0
            push!(non_basic_indices, i)
        end
    end
    return non_basic_indices
end