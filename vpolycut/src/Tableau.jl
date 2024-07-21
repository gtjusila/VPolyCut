#
# Tableau.jl
# A data structure to save the state of the LP tableau 
#
using SCIP

# An LPObject is either an LPColumn or an LPRow
# LPObjects have basis status information

abstract type LPObject end

function get_basis_status(obj::LPObject)::SCIP.SCIP_BASESTAT
    return obj.basis_status
end

# LPColumn Definition and Methods

@kwdef mutable struct LPColumn <: LPObject
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    ub::SCIP.SCIP_Real = -1
    lb::SCIP.SCIP_Real = -1
    sol::SCIP.SCIP_Real = -1
end

function get_col_lb(col::LPColumn)::SCIP.SCIP_Real
    return col.lb
end

function get_col_ub(col::LPColumn)::SCIP.SCIP_Real
    return col.ub
end

function get_col_sol(col::LPColumn)::SCIP.SCIP_Real
    return col.sol
end

#LPRow Definition and Methods

@kwdef mutable struct LPRow <: LPObject
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    rhs::SCIP.SCIP_Real = -1
    lhs::SCIP.SCIP_Real = -1
    slack::SCIP.SCIP_Real = -1
end

function get_row_rhs(row::LPRow)::SCIP.SCIP_Real
    return row.rhs
end

function get_row_lhs(row::LPRow)::SCIP.SCIP_Real
    return row.lhs
end

function get_row_slack(row::LPRow)::SCIP.SCIP_Real
    return row.slack
end

# LPTableau Definition

"""
LPTableau

Concept: SCIP uses a Row Representation of the Tableau. Each LPColumn and LPRow corresponds to a column in the tableau matrix.
Each basic LPColumn or LPRow corresponds to a row in the tableau matrix.

LPTableau inherits the AbstractArray Class and can be indexed as a matrix.
size and getindex are overloaded to access the tableau matrix.
"""
@kwdef mutable struct LPTableau <: Base.AbstractMatrix{SCIP.SCIP_Real}
    # Variables and Rows are distinguished by their unique index given from SCIP
    tableau_matrix::Matrix{SCIP.SCIP_Real} = Matrix{SCIP.SCIP_Real}(undef, 0, 0)
    mapobj2col::Dict{LPObject,Int} = Dict()
    mapcol2obj::Dict{Int,LPObject} = Dict()
    maprow2obj::Dict{Int,LPObject} = Dict()
    mapobj2row::Dict{LPObject,Int} = Dict()
    num_lp_rows::Int = 0
    num_lp_cols::Int = 0
end

# Allow LPTableau to behave as if it is a matrix

function Base.size(A::LPTableau)
    return size(A.tableau_matrix)
end

function Base.getindex(A::LPTableau, i::Int, j::Int)
    return A.tableau_matrix[i, j]
end

function Base.setindex!(A::LPTableau, v, i::Int, j::Int)
    A.tableau_matrix[i, j] = v
end

# Methods for accessing the LPTableau Data

function get_num_tableau_objects(tableau::LPTableau)::Int
    return get_num_lp_cols(tableau) + get_num_lp_rows(tableau)
end

function get_num_basic_objects(tableau::LPTableau)::Int
    return get_num_lp_rows(tableau)
end

function get_num_lp_rows(tableau::LPTableau)::Int
    return tableau.num_lp_rows
end

function get_num_lp_cols(tableau::LPTableau)::Int
    return tableau.num_lp_cols
end

function get_tableau_objects(tableau::LPTableau)::LPObject
    return tableau.objects
end

function get_column_object(tableau::LPTableau, col_index::Int)::LPObject
    return tableau.mapcol2obj[col_index]
end

function get_row_object(tableau::LPTableau, row_index::Int)::LPObject
    return tableau.maprow2obj[row_index]
end

function get_object_column(tableau::LPTableau, obj::LPObject)::Int
    return tableau.mapobj2col[obj]
end

# Methods For Accessing Tableau Data From SCIP

function get_tableau_from_scip(scip::SCIP.SCIPData)::LPTableau
    if SCIP.SCIPisLPSolBasic(scip) != 1
        error("SCIP LP solution is not a basic feasible solution")
    end

    if SCIP.SCIPgetLPSolstat(scip) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        error("SCIP LP solution is not optimal")
    end

    tableau = LPTableau()
    tableau.tableau_matrix = get_full_tableau_matrix(scip)
    tableau.num_lp_cols = SCIP.SCIPgetNLPCols(scip)
    tableau.num_lp_rows = SCIP.SCIPgetNLPRows(scip)

    collect_lp_column_information!(tableau, scip)
    collect_lp_row_information!(tableau, scip)
    collect_basis_infromation!(tableau, scip)

    return tableau
end

function get_full_tableau_matrix(scip::SCIP.SCIPData)::Matrix{SCIP.SCIP_Real}
    # The full tableau is obtained by multiplying the matrix (A I) with B^{-1} 
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    nlprows = SCIP.SCIPgetNLPRows(scip)
    tableau = Matrix{SCIP.SCIP_Real}(undef, nlprows, nlpcols + nlprows)

    for i = 1:nlprows
        tableau[i, :] = get_lp_tableau_row(scip, i)
    end

    return tableau
end

function get_lp_tableau_row(scip::SCIP.SCIPData, row_index::Int)::Vector{SCIP.SCIP_Real}
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    nlprows = SCIP.SCIPgetNLPRows(scip)

    buffer = Vector{SCIP.SCIP_Real}(undef, nlpcols + nlprows)
    # The last `nlprows` entries of the row contains B^{-1}
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvRow(scip, row_index - 1, Ref(buffer, nlpcols + 1), C_NULL, C_NULL)
    # Get the first `nlpcols` entries
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvARow(scip, row_index - 1, C_NULL, Ref(buffer, 1), C_NULL, C_NULL)

    return buffer
end

# Collect LP Columns Information from SCIP

function collect_lp_column_information!(tableau::LPTableau, scip::SCIP.SCIPData)
    lp_columns = get_all_lp_columns(scip)

    for (i, col) in enumerate(lp_columns)
        lp_col = LPColumn()
        lp_col.scip_index = SCIP.SCIPcolGetIndex(col)
        lp_col.basis_status = SCIP.SCIPcolGetBasisStatus(col)
        lp_col.ub = SCIP.SCIPcolGetUb(col)
        lp_col.lb = SCIP.SCIPcolGetLb(col)
        lp_col.sol = SCIP.SCIPcolGetPrimsol(col)

        tableau.mapobj2col[lp_col] = i
        tableau.mapcol2obj[i] = lp_col
    end
end

function get_all_lp_columns(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Col}}
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    lp_columns_ptr = SCIP.SCIPgetLPCols(scip)
    lp_columns = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, lp_columns_ptr, nlpcols)
    return lp_columns
end

# Collect LP Rows Information from SCIP

function collect_lp_row_information!(tableau::LPTableau, scip::SCIP.SCIPData)
    lp_rows = get_all_lp_rows(scip)
    nlpcols = SCIP.SCIPgetNLPCols(scip)

    for (i, row) in enumerate(lp_rows)
        lp_row = LPRow()
        lp_row.scip_index = SCIP.SCIProwGetIndex(row)
        lp_row.basis_status = SCIP.SCIProwGetBasisStatus(row)
        lp_row.rhs = SCIP.SCIProwGetRhs(row)
        lp_row.lhs = SCIP.SCIProwGetLhs(row)
        lp_row.slack = SCIP.SCIPgetRowLPActivity(scip, row)

        col_idx = i + nlpcols
        tableau.mapobj2col[lp_row] = col_idx
        tableau.mapcol2obj[col_idx] = lp_row
    end
end

function get_all_lp_rows(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Row}}
    nlprows = SCIP.SCIPgetNLPRows(scip)
    buffer = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, buffer, nlprows)
    return lp_rows
end

# Collect basis information from SCIP

function collect_basis_infromation!(tableau::LPTableau, scip::SCIP.SCIPData)
    basis_indices = get_basis_indices(scip)
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    for (row, idx) in enumerate(basis_indices)
        if idx < 0
            # Entry is -i, i.e. the ith row is basic
            obj = tableau.mapcol2obj[nlpcols+abs(idx)]
            tableau.maprow2obj[row] = obj
            tableau.mapobj2row[obj] = row
        else
            # Entry is i, i.e. the (i+1)th column is basic (since indexing in C starts with 0)
            obj = tableau.mapcol2obj[idx+1]
            tableau.maprow2obj[row] = obj
            tableau.mapobj2row[obj] = row
        end
    end
end

function get_basis_indices(scip::SCIP.SCIPData)::Vector{Int}
    buffer_length = SCIP.SCIPgetNLPRows(scip)
    buffer = Vector{Cint}(undef, buffer_length)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, Ref(buffer, 1))
    return buffer
end