#
# Tableau.jl
# A data structure to save the state of the LP tableau 
#
using SCIP

"""
Tableau

Suppose we have a linear program in SCIP general form

min c'x
s.t p <= Ax + b <= q
    l <= x <= u

We rewrite the lp into an equivalent problem in the SCIP standard form

min c'x
s.t Ax + b + s = 0
    l <= x <= u
    p <= s <= q

In the SCIP standard form, all of the original constraints have been transformed into 
variables. We can from here on assume that the problem is in SCIP standard form.

The abstract type `Variable` represents every variable in the SCIP standard form. Subsclass
LPRow represents the slack variables s and LPColumn represents the original variables x.

The LP tableau stores tableau information at the optimal solution of the LP in SCIP standard form.
Additionally Tableau may store a constraint matrix that represents the LP in SCIP general form.
"""
@kwdef mutable struct Tableau <: Base.AbstractMatrix{SCIP.SCIP_Real}
    # Variables and Rows are distinguished by their unique index given from SCIP
    tableau_matrix::Matrix{SCIP.SCIP_Real} = Matrix{SCIP.SCIP_Real}(undef, 0, 0)
    mapvar2col::Dict{Tuple{Symbol,Int},Int} = Dict()
    mapcol2var::Dict{Int,Variable} = Dict()
    maprow2var::Dict{Int,Variable} = Dict()
    mapvar2row::Dict{Tuple{Symbol,Int},Int} = Dict()
    noriginalrows::Int = 0
    noriginalcols::Int = 0
    constraint_matrix::Union{Nothing,ConstraintMatrix} = nothing
end

# Allow LPTableau to behave as if it is a matrix
function Base.size(tableau::Tableau)
    return size(tableau.tableau_matrix)
end

function Base.getindex(tableau::Tableau, i::Int, j::Int)
    return tableau.tableau_matrix[i, j]
end

function Base.setindex!(tableau::Tableau, v, i::Int, j::Int)
    tableau.tableau_matrix[i, j] = v
end

# Methods for accessing the LPTableau Data
function has_constraint_matrix_information(tableau::Tableau)::Bool
    return !isnothing(tableau.constraint_matrix)
end

function get_constraint_matrix(tableau::Tableau)::Union{Nothing,ConstraintMatrix}
    if !has_constraint_matrix_information(tableau)
        error("Tableau does not have constraint matrix information")
    end
    return tableau.constraint_matrix
end

function get_nvars(tableau::Tableau)::Int
    return size(tableau)[2]
end

"""
get number of variables in the basis
"""
function get_nbasis(tableau::Tableau)::Int
    return size(tableau)[1]
end

"""
get number of rows in the original problem
"""
function get_noriginalrows(tableau::Tableau)::Int
    return tableau.noriginalrows
end

"""
get number of variables in the original problem
"""
function get_noriginalcols(tableau::Tableau)::Int
    return tableau.noriginalcols
end

function get_var_from_column(tableau::Tableau, col_index::Int)::Variable
    return tableau.mapcol2var[col_index]
end

function get_var_from_row(tableau::Tableau, row_index::Int)::Variable
    return tableau.maprow2var[row_index]
end

function get_column_from_var(tableau::Tableau, var::Variable)::Int
    symbol = get_symbolic_representation(var)
    scip_index = get_scip_index(var)
    return tableau.mapvar2col[(symbol, scip_index)]
end

function set_column_var!(tableau::Tableau, var::Variable, col_index::Int)
    symbol = get_symbolic_representation(var)
    tableau.mapvar2col[(symbol, var.scip_index)] = col_index
    tableau.mapcol2var[col_index] = var
end

function set_row_var!(tableau::Tableau, var::Variable, row_index::Int)
    symbol = get_symbolic_representation(var)
    tableau.mapvar2row[(symbol, var.scip_index)] = row_index
    tableau.maprow2var[row_index] = var
end

function set_constraint_matrix!(tableau::Tableau, matrix::ConstraintMatrix)
    tableau.constraint_matrix = matrix
end

function get_branching_indices(scip::SCIP.SCIPData, tableau::Tableau)::Vector{Int}
    return filter(
        x -> is_branchable(scip, get_var_from_column(tableau, x)), 1:get_nvars(tableau)
    )
end

function is_branchable(scip::SCIP.SCIPData, var::Variable)::Bool
    # If the variable is a slack then we cannot branch
    if isarow(var)
        return false
    end

    # If the variable is not constrainted to be integer we cannot branch
    if SCIP.SCIPvarIsIntegral(get_var_pointer(var)) == 0
        return false
    end

    # If the primal is integral we cannot branch
    if SCIP.SCIPisFeasIntegral(scip, get_sol(var)) == 1
        return false
    end

    return true
end

function convert_standard_inequality_to_general(
    scip::SCIP.SCIPData, tableau::Tableau, standard_row::Vector{SCIP.SCIP_Real}, b
)
    if !has_constraint_matrix_information(tableau)
        error("Tableau does not have constraint matrix information")
    end

    nvars = get_nvars(tableau)
    noriginalcols = get_noriginalcols(tableau)
    general_row = zeros(get_noriginalcols(tableau))
    constraint_matrix = get_constraint_matrix(tableau)

    for i in 1:nvars
        if SCIP.SCIPisZero(scip, standard_row[i]) == 1
            continue
        end
        var = get_var_from_column(tableau, i)
        if isacolumn(var)
            general_row[i] = standard_row[i]
        else
            row_index = i - noriginalcols
            for j in 1:noriginalcols
                general_row[j] -=
                    standard_row[i] * get_entry(constraint_matrix, row_index, j)
            end
            b += standard_row[i] * get_constant(constraint_matrix, row_index)
        end
    end

    return general_row, b
end