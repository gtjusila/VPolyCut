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

const ROW = :ROW
const COLUMN = :COLUMN
# LPColumn Definition and Methods

@kwdef mutable struct LPColumn <: LPObject
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    ub::SCIP.SCIP_Real = -1
    lb::SCIP.SCIP_Real = -1
    sol::SCIP.SCIP_Real = -1
end

function get_lb(col::LPColumn)::SCIP.SCIP_Real
    return col.lb
end

function get_ub(col::LPColumn)::SCIP.SCIP_Real
    return col.ub
end

function get_primal_sol(col::LPColumn)::SCIP.SCIP_Real
    return col.sol
end

#LPRow Definition and Methods
@kwdef mutable struct LPRow <: LPObject
    scip_index::Int = -1
    basis_status::SCIP.SCIP_BASESTAT = SCIP.SCIP_BASESTAT_ZERO
    rhs::SCIP.SCIP_Real = -1
    lhs::SCIP.SCIP_Real = -1
    obj_pointer::Ptr{SCIP.SCIP_ROW} = C_NULL
    slack::SCIP.SCIP_Real = -1
end

function get_rhs(row::LPRow)::SCIP.SCIP_Real
    return row.rhs
end

function get_lhs(row::LPRow)::SCIP.SCIP_Real
    return row.lhs
end

function get_slack(row::LPRow)::SCIP.SCIP_Real
    return row.slack
end

function get_pointer(row::LPRow)::Ptr{SCIP.SCIP_ROW}
    return row.obj_pointer
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
    mapobj2col::Dict{Tuple{Symbol,Int},Int} = Dict()
    mapcol2obj::Dict{Int,LPObject} = Dict()
    maprow2obj::Dict{Int,LPObject} = Dict()
    mapobj2row::Dict{Tuple{Symbol,Int},Int} = Dict()
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

function get_nobjects(tableau::LPTableau)::Int
    return get_nlpcols(tableau) + get_nlprows(tableau)
end

function get_nbasicobjects(tableau::LPTableau)::Int
    return get_nlprows(tableau)
end

function get_nlprows(tableau::LPTableau)::Int
    return tableau.num_lp_rows
end

function get_nlpcols(tableau::LPTableau)::Int
    return tableau.num_lp_cols
end

function get_objects(tableau::LPTableau)::LPObject
    return tableau.objects
end

function get_column_object(tableau::LPTableau, col_index::Int)::LPObject
    return tableau.mapcol2obj[col_index]
end

function get_row_object(tableau::LPTableau, row_index::Int)::LPObject
    return tableau.maprow2obj[row_index]
end

function get_column_index(tableau::LPTableau, obj::LPObject)::Int
    if isa(obj, LPColumn)
        return tableau.mapobj2col[(COLUMN, obj.scip_index)]
    else
        return tableau.mapobj2col[(ROW, obj.scip_index)]
    end
end

function get_object_from_scip_index(tableau::LPTableau, scip_index::Int, object_type::Symbol)::LPObject
    # Create a dummy object and use it to look for the actual object within the tableau columns
    dummy = nothing
    if object_type == COLUMN
        dummy = LPColumn(scip_index, SCIP.SCIP_BASESTAT_ZERO, -1, -1, -1)
    else
        dummy = LPRow(scip_index, SCIP.SCIP_BASESTAT_ZERO, -1, -1, C_NULL, -1)
    end
    col_index = get_column_index(tableau, dummy)

    return get_column_object(tableau, col_index)
end

# Setter Methods for tableau
function set_column_object!(tableau::LPTableau, obj::LPObject, col_index::Int)
    if isa(obj, LPColumn)
        tableau.mapobj2col[(COLUMN, obj.scip_index)] = col_index
    elseif isa(obj, LPRow)
        tableau.mapobj2col[(ROW, obj.scip_index)] = col_index
    end
    tableau.mapcol2obj[col_index] = obj
end

function set_row_object!(tableau::LPTableau, obj::LPObject, row_index::Int)
    if isa(obj, LPColumn)
        tableau.mapobj2row[(COLUMN, obj.scip_index)] = row_index
    elseif isa(obj, LPRow)
        tableau.mapobj2row[(ROW, obj.scip_index)] = row_index
    end
    tableau.maprow2obj[row_index] = obj
end

# Methods For Accessing Tableau Data From SCIP

function construct_tableau(scip::SCIP.SCIPData)::LPTableau
    if SCIP.SCIPisLPSolBasic(scip) != 1
        error("SCIP LP solution is not a basic feasible solution")
    end

    if SCIP.SCIPgetLPSolstat(scip) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        error("SCIP LP solution is not optimal")
    end

    tableau = LPTableau()
    tableau.tableau_matrix = fetch_tableau_matrix(scip)
    tableau.num_lp_cols = SCIP.SCIPgetNLPCols(scip)
    tableau.num_lp_rows = SCIP.SCIPgetNLPRows(scip)

    fetch_columns!(tableau, scip)
    fetch_rows!(tableau, scip)
    fetch_basis_information!(tableau, scip)

    return tableau
end

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

# Collect LP Columns Information from SCIP

function fetch_columns!(tableau::LPTableau, scip::SCIP.SCIPData)
    lp_columns = fetch_column_pointers(scip)

    for (i, col) in enumerate(lp_columns)
        lp_col = LPColumn()
        lp_col.scip_index = SCIP.SCIPcolGetIndex(col)
        lp_col.basis_status = SCIP.SCIPcolGetBasisStatus(col)
        lp_col.ub = SCIP.SCIPcolGetUb(col)
        lp_col.lb = SCIP.SCIPcolGetLb(col)
        lp_col.sol = SCIP.SCIPcolGetPrimsol(col)

        set_column_object!(tableau, lp_col, i)
    end
end

function fetch_column_pointers(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Col}}
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    lp_columns_ptr = SCIP.SCIPgetLPCols(scip)
    lp_columns = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Col}}, lp_columns_ptr, nlpcols)
    return lp_columns
end

# Collect LP Rows Information from SCIP

function fetch_rows!(tableau::LPTableau, scip::SCIP.SCIPData)
    lp_rows = fetch_row_pointers(scip)
    nlpcols = SCIP.SCIPgetNLPCols(scip)

    for (i, row) in enumerate(lp_rows)
        lp_row = LPRow()
        lp_row.scip_index = SCIP.SCIProwGetIndex(row)
        lp_row.basis_status = SCIP.SCIProwGetBasisStatus(row)
        lp_row.rhs = SCIP.SCIProwGetRhs(row)
        lp_row.lhs = SCIP.SCIProwGetLhs(row)
        lp_row.slack = SCIP.SCIPgetRowLPActivity(scip, row)
        lp_row.obj_pointer = row

        if SCIP.SCIProwGetConstant(row) != 0
            @warn "Row has non-zero constant term"
        end

        col_idx = i + nlpcols
        set_column_object!(tableau, lp_row, col_idx)
    end
end

function fetch_row_pointers(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_Row}}
    nlprows = SCIP.SCIPgetNLPRows(scip)
    buffer = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_Row}}, buffer, nlprows)
    return lp_rows
end

# Collect basis information from SCIP

function fetch_basis_information!(tableau::LPTableau, scip::SCIP.SCIPData)
    basis_indices = fetch_basic_indices(scip)
    nlpcols = SCIP.SCIPgetNLPCols(scip)
    for (row, idx) in enumerate(basis_indices)
        if idx < 0
            # Entry is -i, i.e. the ith row is basic
            obj = tableau.mapcol2obj[nlpcols+abs(idx)]
        else
            # Entry is i, i.e. the (i+1)th column is basic (since indexing in C starts with 0)
            obj = tableau.mapcol2obj[idx+1]
        end
        set_row_object!(tableau, obj, row)
    end
end

function fetch_basic_indices(scip::SCIP.SCIPData)::Vector{Int}
    buffer_length = SCIP.SCIPgetNLPRows(scip)
    buffer = Vector{Cint}(undef, buffer_length)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, Ref(buffer, 1))
    return buffer
end