#
# Tableau.jl
# A data structure to save the state of the LP tableau 
#
using SCIP

"""
LPTableau

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
function get_nvars(tableau::Tableau)::Int
    return size(tableau.tableau_matrix)[2]
end

"""
get number of variables in the basis
"""
function get_nbasis(tableau::Tableau)::Int
    return size(tableau.tableau_matrix)[1]
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

function get_all_variables(tableau::Tableau)::Variable
    return collect(values(tableau.mapcol2var))
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

#TODO CLEANUP
function convert_standard_row_to_general(tableau::Tableau, standard_row::Vector{SCIP.SCIP_Real}, b)
    general_row = zeros(get_noriginalcols(tableau))
    for i = 1:get_nvars(tableau)
        if standard_row[i] != 0
            if i > get_noriginalcols(tableau)
                original_row = get_var_from_column(tableau, i)
                original_row_coefficients = get_row_coefficients(original_row)
                b += standard_row[i] * get_constant(original_row)
                for j = 1:get_noriginalcols(tableau)
                    general_row[j] -= standard_row[i] * original_row_coefficients[j]
                end
            else
                general_row[i] = standard_row[i]
            end
        end
    end
    return general_row, b
end