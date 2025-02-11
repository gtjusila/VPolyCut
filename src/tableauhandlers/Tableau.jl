#
# Tableau.jl
# A data structure to save the state of the LP tableau 
#
using SCIP

"""
Tableau

The tableau class allows the user to temporarily "freeze" the solution of the SCIP LP interface.
Suppose we have a linear program in SCIP general form

min c'x
s.t p <= Ax + b <= q
    l <= x <= u

We rewrite the lp into an equivalent problem in the SCIP standard form

min c'x
s.t Ax + b - s = 0
    l <= x <= u
    p <= s <= q

In the SCIP standard form, all of the original inequality constraints have been transformed into 
equality constraint with slack variables. The motivation behind defining the SCIP standard form this way is because
SCIP_BASESTAT_UPPER for a row means that s=q and SCIP_BASESTAT_LOWER for a row
means that s=p (that is the row is tight at the upper and lower bound respectively). 
There is however a twist. SCIP LPI requires that the coefficient of slack variables s has to be +1. 
Defining r = -s, we have

min c'x
s.t. Ax + b + r = 0
     l <= x <= u
     -q <= r <= -p

This formulation satisfies the formulation requirement of the LPI. We call this form the SCIP 
the SCIP Normal Form. The struct Tableau assumes the problem is in SCIP Normal Form

The abstract type `Variable` represents every variable in the SCIP normal form. Subsclass
LPRow represents the normalized slack variables r and LPColumn represents the original variables x.

The LP tableau stores tableau information at the optimal solution of the LP in SCIP normal form.
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

"""
get_nvars

Get the number of variables in the tableau (original and slack variables).
This is equal to the sum of the number of inequality and the number of variables in the original problem.
"""
get_nvars(tableau::Tableau) = size(tableau)[2]

"""
get_nbasis

Get the number of variables in the basis (i.e. number of basic variables). This is equal to the number of inequality constraints 
"""
get_nbasis(tableau::Tableau) = size(tableau)[1] # Get number of variables in the basis

"""
get_noriginalrows

Get the number of rows in the original problem (i.e. number of inequality constraints)
"""
get_noriginalrows(tableau::Tableau) = tableau.noriginalrows # Get number of rows in the original problem

"""
get_noriginalcols

Get the number of columns in the original problem (i.e. number of variables in the original problem)
"""
get_noriginalcols(tableau::Tableau) = tableau.noriginalcols # Get number of variables in the original problem

"""
Base.size

The size of a tableau can also be prompted using the julia size function. The size of the tableau matrix is
(number of inequalities) * (number of variables + number of inequalites)
"""
Base.size(tableau::Tableau) = size(tableau.tableau_matrix)

"""
Base.getindex

get the value of the tableau matrix at position (i, j)
"""
Base.getindex(tableau::Tableau, i::Int, j::Int) = tableau.tableau_matrix[i, j]

"""
Base.setindex!
"""
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
    if is_integral(get_sol(var))
        return false
    end

    return true
end

function convert_standard_inequality_to_general(
    tableau::Tableau, standard_row::Vector{SCIP.SCIP_Real}, b
)
    if !has_constraint_matrix_information(tableau)
        error("Tableau does not have constraint matrix information")
    end

    nvars = get_nvars(tableau)
    noriginalcols = get_noriginalcols(tableau)
    general_row = zeros(get_noriginalcols(tableau))
    constraint_matrix = get_constraint_matrix(tableau)

    for i in 1:nvars
        if is_zero(standard_row[i])
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

"""
The Objective function is stored such that it is a function of nonbasic variables 
"""
function get_objective_direction(tableau::Tableau)::Vector{SCIP.SCIP_Real}
    dim = get_nvars(tableau)
    objective_direction = zeros(dim)
    for i in 1:dim
        var = get_var_from_column(tableau, i)
        objective_direction[i] = get_obj(var)
    end
    return objective_direction
end

"""
Get Pointers to problem variables ordered by their column index in the tableau
"""
function get_problem_variables_pointers(tableau::Tableau)::Vector{Ptr{SCIP.SCIP_VAR}}
    noriginalcols = get_noriginalcols(tableau)
    pointers = Vector{Ptr{SCIP.SCIP_VAR}}(undef, noriginalcols)
    for i in 1:noriginalcols
        pointers[i] = get_var_pointer(get_var_from_column(tableau, i))
    end
    return pointers
end

function get_tableau_density(scip::SCIP.SCIPData, tableau::Tableau)::Float64
    return sum([!is_zero(i) ? 1 : 0 for i in tableau.tableau_matrix]) /
           length(tableau.tableau_matrix)
end