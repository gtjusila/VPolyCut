# This file contains the logic for extracting tableau data from SCIP 

# Methods For Accessing Tableau Data From SCIP
function construct_tableau(scip::SCIP.SCIPData)::Tableau
    if SCIP.SCIPisLPSolBasic(scip) != 1
        error("SCIP LP solution is not a basic feasible solution")
    end

    if SCIP.SCIPgetLPSolstat(scip) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        error("SCIP LP solution is not optimal")
    end

    tableau = Tableau()
    tableau.tableau_matrix = fetch_tableau_matrix(scip)
    tableau.noriginalcols = SCIP.SCIPgetNLPCols(scip)
    tableau.noriginalrows = SCIP.SCIPgetNLPRows(scip)

    fetch_columns!(tableau, scip)
    fetch_rows!(tableau, scip)
    fetch_basis_information!(tableau, scip)

    return tableau
end

# Code for fetching tableau matrix from SCIP

function fetch_tableau_matrix(scip::SCIP.SCIPData)::Matrix{SCIP.SCIP_Real}
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
function fetch_columns!(tableau::Tableau, scip::SCIP.SCIPData)
    lp_columns = fetch_column_pointers(scip)

    for (i, col) in enumerate(lp_columns)
        lp_col = LPColumn()
        populate!(scip, col, lp_col)
        set_column_var!(tableau, lp_col, i)
    end
end

function fetch_column_pointers(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Col}}
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    lp_columns_ptr = SCIP.SCIPgetLPCols(scip)
    lp_columns = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, lp_columns_ptr, nlpcols)
    return lp_columns
end

# Collect LP Rows Information from SCIP
#TODO CLEANUP
function fetch_rows!(tableau::Tableau, scip::SCIP.SCIPData)
    lp_rows = fetch_row_pointers(scip)
    nlpcols = SCIP.SCIPgetNLPCols(scip)

    # Fetch Rows must be done after fetch columns
    @assert length(tableau.mapcol2var) == nlpcols
    for (i, row) in enumerate(lp_rows)
        lp_row = LPRow()
        populate!(scip, row, lp_row)

        coef = get_row_coefficients(row, tableau)
        set_row_coefficient!(lp_row, coef)
        if SCIP.SCIProwGetConstant(row) != 0
            @warn "Row has non-zero constant term"
        end

        col_idx = i + nlpcols
        set_column_var!(tableau, lp_row, col_idx)
    end
end

# TODO CLEANUP
function get_row_coefficients(scip_ptr::Ptr{SCIP.SCIP_ROW}, tableau::Tableau)::Vector{SCIP.SCIP_Real}
    # Get the coefficients of the row
    nnonz = SCIP.SCIProwGetNNonz(scip_ptr)
    vars_ptr = SCIP.SCIProwGetCols(scip_ptr)
    vars_ptr = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, vars_ptr, nnonz)
    vars = Variable[]
    for v in vars_ptr
        # Here we unfortunately need to do low level manipulation to get the Variable
        col_symbol = get_symbolic_representation(LPColumn())
        col = tableau.mapvar2col[(col_symbol, Int(SCIP.SCIPcolGetIndex(v)))]
        var = tableau.mapcol2var[col]
        push!(vars, var)
    end
    coefs = SCIP.SCIProwGetVals(scip_ptr)
    coefs = unsafe_wrap(Vector{SCIP.SCIP_Real}, coefs, nnonz)
    dim = get_noriginalcols(tableau)
    row = zeros(dim)
    for i = 1:nnonz
        row[get_column_from_var(tableau, vars[i])] = coefs[i]
    end
    return row
end

function fetch_row_pointers(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Row}}
    nlprows = SCIP.SCIPgetNLPRows(scip)
    buffer = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, buffer, nlprows)
    return lp_rows
end

# Collect basis information from SCIP

function fetch_basis_information!(tableau::Tableau, scip::SCIP.SCIPData)
    basis_indices = fetch_basic_indices(scip)
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

function fetch_basic_indices(scip::SCIP.SCIPData)::Vector{Int}
    buffer_length = SCIP.SCIPgetNLPRows(scip)
    buffer = Vector{Cint}(undef, buffer_length)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, Ref(buffer, 1))
    return buffer
end