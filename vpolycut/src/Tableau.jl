using SCIP

@kwdef mutable struct TableauColumn
    index::Int = -1
    scip_pointer::Ptr{SCIP.SCIP_Col} = C_NULL
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    upper_bound::SCIP.SCIP_Real = -1
    lower_bound::SCIP.SCIP_Real = -1
    tableau_column::Int = -1
    tableau_row::Int = -1
end

@kwdef mutable struct TableauRow
    index::Int = -1
    scip_pointer::Ptr{SCIP.SCIP_Row} = C_NULL
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    rhs::SCIP.SCIP_Real = -1
    lhs::SCIP.SCIP_Real = -1
    tableau_column::Int = -1
    tableau_row::Int = -1
end

"""
DenseTableau

Concept: TableauColumns and TableauRows represent variables in the LP problem. 
TableauColumns correspond to the problem variables, while TableauRows correspond to the slack variables for the constraints. 
Each column of the tableau matrix corresponds to either an TableauColumn or an TableauRow. Similarly, 
each row in the matrix corresponds to one basic TableauColumn or TableauRow.

The column index in the tableau matrix for an TableauColumn or TableauRow can be obtained using:
- get_matrix_column_index_for_tableau_column
- get_matrix_column_index_for_tableu_row

If the TableauColumn or TableauRow is basic, the corresponding row index in the tableau matrix can be obtained using:
- get_matrix_row_index_for_tableau_column
- get_matrix_row_index_for_tableau_row
"""

@kwdef mutable struct DenseTableau
    # Variables and Rows are distinguished by their unique index given from SCIP
    tableau_matrix::Matrix{SCIP.SCIP_Real} = Matrix{SCIP.SCIP_Real}(undef, 1, 1)
    lp_columns::Dict{Int,TableauColumn} = Dict()
    lp_rows::Dict{Int,TableauRow} = Dict()
end

#
# Methods to access tableau data
#

function get_tableau_columns(tableau::DenseTableau)::Vector{TableauColumn}
    return collect(values(tableau.lp_columns))
end

function get_tableau_column_basis_status(col::TableauColumn)::SCIP.SCIP_BASESTAT
    return col.basis_status
end

function get_tableau_column_upper_bound(col::TableauColumn)::SCIP.SCIP_Real
    return col.upper_bound
end

function get_tableau_column_lower_bound(col::TableauColumn)::SCIP.SCIP_Real
    return col.lower_bound
end

function get_matrix_column_index_for_tableau_column(col::TableauColumn)::Int
    return col.tableau_column
end

function get_matrix_row_index_for_tableau_column(col::TableauColumn)::Int
    if get_tableau_column_basis_status(col) == SCIP.SCIP_BASESTAT_BASIC
        return col.tableau_row
    else
        return -1
    end
end

function get_tableau_rows(tableau::DenseTableau)::Vector{TableauRow}
    return collect(values(tableau.lp_rows))
end

function get_tableau_row_basis_status(row::TableauRow)::SCIP.SCIP_BASESTAT
    return row.basis_status
end

function get_tableau_row_rhs(row::TableauRow)::SCIP.SCIP_Real
    return row.rhs
end

function get_tableau_row_lhs(row::TableauRow)::SCIP.SCIP_Real
    return row.lhs
end

function get_matrix_column_index_for_tableau_row(row::TableauRow)::Int
    return row.tableau_column
end

function get_matrix_row_index_for_tableau_row(row::TableauRow)::Int
    if get_tableau_row_basis_status(row) == SCIP.SCIP_BASESTAT_BASIC
        return row.tableau_row
    else
        return -1
    end
end

#
# Methods creating and filling in DenseTableau 
#
function get_dense_tableau_from_scip(scip::SCIP.SCIPData)::DenseTableau
    @assert SCIP.SCIPisLPSolBasic(scip) == 1
    @assert SCIP.SCIPgetLPSolstat(scip) == SCIP.SCIP_LPSOLSTAT_OPTIMAL

    tableau = DenseTableau()
    tableau.tableau_matrix = get_full_tableau_matrix(scip)

    collect_column_information!(tableau, scip)
    collect_row_information!(tableau, scip)

    return tableau
end

function get_full_tableau_matrix(scip::SCIP.SCIPData)::Matrix{SCIP.SCIP_Real}
    # The full tableau is obtained by multiplying the matrix (A I) with B^{-1} 
    var_count = get_lp_column_count(scip)
    row_count = get_lp_row_count(scip)
    tableau = Matrix{SCIP.SCIP_Real}(undef, row_count, var_count + row_count)

    for i = 1:row_count
        tableau[i, :] = get_lp_tableau_row(scip, i)
    end

    return tableau
end

function collect_column_information!(tableau::DenseTableau, scip::SCIP.SCIPData)
    lp_columns = get_all_lp_columns(scip)

    for col in lp_columns
        tableau_col = TableauColumn()
        col_index = get_lp_column_index(col)

        tableau_col.index = col_index
        tableau_col.scip_pointer = col
        tableau_col.basis_status = get_lp_column_basis_status(col)
        tableau_col.upper_bound = get_lp_column_upper_bound(col)
        tableau_col.lower_bound = get_lp_column_lower_bound(col)
        tableau_col.tableau_column = get_lp_column_position(col)

        if is_lp_column_basic(col)
            tableau_col.tableau_row = get_lp_column_position_in_basis(scip, col)
        end

        tableau.lp_columns[col_index] = tableau_col
    end
end

function get_basis_indices(scip::SCIP.SCIPData)::Vector{Int}
    buffer_length = get_lp_row_count(scip)
    buffer = Vector{Cint}(undef, buffer_length)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, Ref(buffer, 1))
    return buffer
end

function get_lp_column_upper_bound(col::Ptr{SCIP.SCIP_Col})::SCIP.SCIP_Real
    return SCIP.SCIPcolGetUb(col)
end

function get_lp_column_lower_bound(col::Ptr{SCIP.SCIP_Col})::SCIP.SCIP_Real
    return SCIP.SCIPcolGetLb(col)
end

function get_lp_column_position_in_basis(scip::SCIP.SCIPData, col::Ptr{SCIP.SCIP_Col})::Int
    @assert is_lp_column_basic(col)
    basis_indices = get_basis_indices(scip)
    col_pos = get_lp_column_position(col)
    return findfirst(isequal(col_pos - 1), basis_indices)
end

function is_lp_column_basic(col::Ptr{SCIP.SCIP_Col})::Bool
    return (get_lp_column_basis_status(col) == SCIP.SCIP_BASESTAT_BASIC)
end

function get_lp_column_basis_status(col::Ptr{SCIP.SCIP_Col})::SCIP.SCIP_BASESTAT
    return SCIP.SCIPcolGetBasisStatus(col)
end

function get_lp_column_position(col::Ptr{SCIP.SCIP_Col})::Int
    return SCIP.SCIPcolGetLPPos(col) + 1
end

function get_lp_column_count(scip::SCIP.SCIPData)::Int
    return SCIP.SCIPgetNLPCols(scip)
end

function get_lp_column_index(col::Ptr{SCIP.SCIP_Col})::Int
    return SCIP.SCIPcolGetIndex(col)
end

function get_all_lp_columns(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Col}}
    column_count = get_lp_column_count(scip)
    lp_columns_ptr = SCIP.SCIPgetLPCols(scip)
    lp_columns = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, lp_columns_ptr, column_count)
    return lp_columns
end

function collect_row_information!(tableau::DenseTableau, scip::SCIP.SCIPData)
    lp_rows = get_all_lp_rows(scip)
    for row in lp_rows
        tableau_row = TableauRow()
        row_index = get_lp_row_index(row)

        tableau_row.index = row_index
        tableau_row.scip_pointer = row
        tableau_row.rhs = get_lp_row_rhs(row)
        tableau_row.lhs = get_lp_row_lhs(row)
        tableau_row.basis_status = get_lp_row_basis_status(row)
        tableau_row.tableau_column = get_lp_column_count(scip) + get_lp_row_position(row)

        if is_lp_row_basic(row)
            tableau_row.tableau_row = get_lp_row_position_in_basis(scip, row)
        end

        tableau.lp_rows[row_index] = tableau_row
    end
end

function get_lp_row_rhs(row::Ptr{SCIP.SCIP_Row})::SCIP.SCIP_Real
    return SCIP.SCIProwGetRhs(row)
end

function get_lp_row_lhs(row::Ptr{SCIP.SCIP_Row})::SCIP.SCIP_Real
    return SCIP.SCIProwGetLhs(row)
end

function get_lp_row_position_in_basis(scip::SCIP.SCIPData, row::Ptr{SCIP.SCIP_Row})::Int
    @assert is_lp_row_basic(row)
    basis_indices = get_basis_indices(scip)
    row_pos = get_lp_row_position(row)
    return findfirst(isequal(-1 * row_pos), basis_indices)
end

function is_lp_row_basic(row::Ptr{SCIP.SCIP_Row})::Bool
    return (get_lp_row_basis_status(row) == SCIP.SCIP_BASESTAT_BASIC)
end

function get_lp_row_basis_status(row::Ptr{SCIP.SCIP_Row})::SCIP.SCIP_BASESTAT
    return SCIP.SCIProwGetBasisStatus(row)
end

function get_lp_row_position(row::Ptr{SCIP.SCIP_Row})::Int
    return SCIP.SCIProwGetLPPos(row) + 1
end

function get_lp_row_count(scip::SCIP.SCIPData)::Int
    return SCIP.SCIPgetNLPRows(scip)
end

function get_lp_row_index(row::Ptr{SCIP.SCIP_Row})::Int
    return SCIP.SCIProwGetIndex(row)
end

function get_all_lp_rows(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Row}}
    row_count = get_lp_row_count(scip)
    buffer = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, buffer, row_count)
    return lp_rows
end

function get_lp_tableau_row(scip::SCIP.SCIPData, row_index::Int)::Vector{SCIP.SCIP_Real}
    column_count = get_lp_column_count(scip)
    row_count = get_lp_row_count(scip)
    buffer = Vector{SCIP.SCIP_Real}(undef, column_count + row_count)

    # The last `row_count` entries of the row contains B^{-1}
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvRow(scip, row_index - 1, Ref(buffer, column_count + 1), C_NULL, C_NULL)

    # Get the first `var_count` entries
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvARow(scip, row_index - 1, C_NULL, Ref(buffer, 1), C_NULL, C_NULL)

    return buffer
end

