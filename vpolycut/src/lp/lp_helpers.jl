using SCIP

#
function get_scip_tableau(scip::SCIP.SCIPData)
    row_count, column_count = get_scip_tableau_dimension(scip)
    tableau = create_tableau(row_count, column_count)

    for i = 1:row_count
        # Allocate Memory to store tableau entries
        row_vector = Vector{SCIP.SCIP_Real}(undef, column_count)
        get_scip_tableau_ith_row(scip, i, row_vector)
        copy_vector_to_tableau_row(tableau, i, row_vector)
    end

    return tableau
end

function get_scip_tableau_ith_row(scip::SCIP.SCIPData, i::UInt, buffer::Vector{SCIP.SCIP_Real})
    # Warning! C-indexing/Julia index mix coming
    # Tableau is given by [B^{-1}A B^{-1}]
    var_count = SCIP.SCIPgetNLPCols(scip)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvRow(scip, i - 1, Ref(buffer, var_count + 1), C_NULL, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvARow(scip, i - 1, Ref(buffer, var_count + 1), Ref(buffer, 1), C_NULL, C_NULL)
end

function get_scip_tableau_dimension(scip::SCIP.SCIPData)::Tuple{UInt,UInt}
    row_count = SCIP.SCIPgetNLPRows(scip)
    column_count = SCIP.SCIPgetNLPCols(scip) # LP have N Columns
    return (row_count, row_count + column_count)
end
