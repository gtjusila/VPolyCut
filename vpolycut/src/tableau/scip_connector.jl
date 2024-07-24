# This file contains the logic for extracting tableau data from SCIP 

# Methods For Accessing Tableau Data From SCIP
function construct_tableau_with_constraint_matrix(scip::SCIP.SCIPData)::Tableau
    tableau = construct_tableau(scip)
    constraint_matrix = fetch_constraint_matrix_from_scip(tableau, scip)
    set_constraint_matrix!(tableau, constraint_matrix)
    return tableau
end

function construct_tableau(scip::SCIP.SCIPData)::Tableau
    if SCIP.SCIPisLPSolBasic(scip) != 1
        error("SCIP LP solution is not a basic feasible solution")
    end

    if SCIP.SCIPgetLPSolstat(scip) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        error("SCIP LP solution is not optimal")
    end

    tableau = Tableau()
    tableau.tableau_matrix = fetch_tableau_matrix_from_scip(scip)
    tableau.noriginalcols = SCIP.SCIPgetNLPCols(scip)
    tableau.noriginalrows = SCIP.SCIPgetNLPRows(scip)

    add_variables_from_scip_columns!(tableau, scip)
    add_variables_from_scip_rows!(tableau, scip)
    fetch_basis_data_from_scip!(tableau, scip)

    return tableau
end

# Code for fetching tableau matrix from SCIP
function fetch_tableau_matrix_from_scip(scip::SCIP.SCIPData)::Matrix{SCIP.SCIP_Real}
    # The full tableau is obtained by multiplying the matrix (A I) with B^{-1} 
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    nlprows = SCIP.SCIPgetNLPRows(scip)
    tableau = Matrix{SCIP.SCIP_Real}(undef, nlprows, nlpcols + nlprows)

    for i = 1:nlprows
        tableau[i, :] = fetch_tableau_row(scip, i)
    end

    return tableau
end

function fetch_tableau_row(scip::SCIP.SCIPData, row_index::Int)::Vector{SCIP.SCIP_Real}
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    nlprows = SCIP.SCIPgetNLPRows(scip)

    buffer = Vector{SCIP.SCIP_Real}(undef, nlpcols + nlprows)
    # The last `nlprows` entries of the row contains B^{-1}
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvRow(scip, row_index - 1, Ref(buffer, nlpcols + 1), C_NULL, C_NULL)
    # Get the first `nlpcols` entries
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvARow(scip, row_index - 1, C_NULL, Ref(buffer, 1), C_NULL, C_NULL)

    return buffer
end

#Collect LP Columns Information from SCIP
function add_variables_from_scip_columns!(tableau::Tableau, scip::SCIP.SCIPData)
    cols_pointers = fetch_column_pointers(scip)

    for col_pointer in cols_pointers
        var = create_variable_from_column_pointer(col_pointer)
        # Casted because SCIP returns Int32 for indices
        # Plus 1 because C index start with 0
        index::Int = SCIP.SCIPcolGetLPPos(col_pointer) + 1
        set_column_var!(tableau, var, index)
    end
end

function create_variable_from_column_pointer(scip_ptr::Ptr{SCIP.SCIP_COL})::Variable
    lp_col = LPColumn()
    set_scip_index!(lp_col, SCIP.SCIPcolGetIndex(scip_ptr))
    set_basis_status!(lp_col, SCIP.SCIPcolGetBasisStatus(scip_ptr))
    set_ub!(lp_col, SCIP.SCIPcolGetUb(scip_ptr))
    set_lb!(lp_col, SCIP.SCIPcolGetLb(scip_ptr))
    set_sol!(lp_col, SCIP.SCIPcolGetPrimsol(scip_ptr))
    set_var_pointer!(lp_col, SCIP.SCIPcolGetVar(scip_ptr))
    return lp_col
end

function fetch_column_pointers(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Col}}
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    lp_cols_ptr = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, lp_cols_ptr, nlpcols)
    return lp_cols
end

# Collect LP Rows Information from SCIP
function add_variables_from_scip_rows!(tableau::Tableau, scip::SCIP.SCIPData)
    row_pointers = fetch_row_pointers(scip)
    nlpcols = SCIP.SCIPgetNLPCols(scip)

    for row_pointer in row_pointers
        var = create_variable_from_row_pointer(scip, row_pointer)
        # Column variables come before row variables so we offset by nlpcols
        # Casted because SCIP returns Int32 for indices
        # Plus 1 because C index start with 0
        index::Int = nlpcols + SCIP.SCIProwGetLPPos(row_pointer) + 1
        set_column_var!(tableau, var, index)
    end
end

function create_variable_from_row_pointer(scip::SCIP.SCIPData, scip_ptr::Ptr{SCIP.SCIP_ROW})::Variable
    lp_row = LPRow()
    set_scip_index!(lp_row, SCIP.SCIProwGetIndex(scip_ptr))
    set_basis_status!(lp_row, SCIP.SCIProwGetBasisStatus(scip_ptr))
    set_ub!(lp_row, SCIP.SCIProwGetRhs(scip_ptr))
    set_lb!(lp_row, SCIP.SCIProwGetLhs(scip_ptr))
    set_sol!(lp_row, -SCIP.SCIPgetRowActivity(scip, scip_ptr))
    return lp_row
end

function fetch_row_pointers(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Row}}
    nlprows = SCIP.SCIPgetNLPRows(scip)
    buffer = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, buffer, nlprows)
    return lp_rows
end

# Collect basis information from SCIP
function fetch_basis_data_from_scip!(tableau::Tableau, scip::SCIP.SCIPData)
    basis_indices = fetch_basis_indices(scip)
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    for (row, idx) in enumerate(basis_indices)
        if idx < 0
            # Entry is -i, i.e. the ith row is basic
            var = tableau.mapcol2var[nlpcols+abs(idx)]
        else
            # Entry is i, i.e. the (i+1)th column is basic (since indexing in C starts with 0)
            var = tableau.mapcol2var[idx+1]
        end
        set_row_var!(tableau, var, row)
    end
end

function fetch_basis_indices(scip::SCIP.SCIPData)::Vector{Int}
    buffer_length = SCIP.SCIPgetNLPRows(scip)
    buffer = Vector{Cint}(undef, buffer_length)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, Ref(buffer, 1))
    return buffer
end

function fetch_constraint_matrix_from_scip(tableau::Tableau, scip::SCIP.SCIPData)::ConstraintMatrix
    # Get the coefficients of the row
    noriginalrows = get_noriginalrows(tableau)
    noriginalcols = get_noriginalcols(tableau)
    matrix = ConstraintMatrix(noriginalrows, noriginalcols)

    row_pointers = fetch_row_pointers(scip)
    for (row_index, row_pointer) in enumerate(row_pointers)
        row_entries = fetch_non_zero_entries_from_row(row_pointer)
        for (col_index, entry) in row_entries
            set_entry!(matrix, row_index, col_index, entry)
        end
        set_constant!(matrix, row_index, SCIP.SCIProwGetConstant(row_pointer))
    end

    return matrix
end

function fetch_non_zero_entries_from_row(row_pointer::Ptr{SCIP.SCIP_ROW})::Vector{Tuple{Int,SCIP.SCIP_Real}}
    buffer::Vector{Tuple{Int,SCIP.SCIP_Real}} = []
    nnonzeros = SCIP.SCIProwGetNNonz(row_pointer)

    nonzero_columns = SCIP.SCIProwGetCols(row_pointer)
    nonzero_columns = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, nonzero_columns, nnonzeros)
    nonzero_entries = SCIP.SCIProwGetVals(row_pointer)
    nonzero_entries = unsafe_wrap(Vector{SCIP.SCIP_Real}, nonzero_entries, nnonzeros)

    for i = 1:nnonzeros
        col_index = SCIP.SCIPcolGetLPPos(nonzero_columns[i]) + 1
        col_entry = nonzero_entries[i]
        if col_index > 0
            push!(buffer, (col_index, col_entry))
        end
    end

    return buffer
end



