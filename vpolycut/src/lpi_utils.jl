import SCIP
include("constants.jl")

"""
Take a vector alpha of the problem dimension and add row lhs <= alpha*x <= rhs to LPI
"""
function add_row_to_lpi(
    lpi::Ptr{SCIP.SCIP_LPi}; 
    alpha::Vector{SCIP.SCIP_Real},
    lhs::Union{SCIP.SCIP_Real, Nothing} = nothing,
    rhs::Union{SCIP.SCIP_Real, Nothing} = nothing,
    row_name::String = ""
)
    @assert(!(isnothing(lhs) && isnothing(rhs))) #LHS and RHS cannot both be nothing
    if isnothing(lhs) lhs = -SCIP.SCIPlpiInfinity(lpi) end
    if isnothing(rhs) rhs = SCIP.SCIPlpiInfinity(lpi) end
   
    lhs = Ref{SCIP.SCIP_Real}(lhs)
    rhs = Ref{SCIP.SCIP_Real}(rhs)
    
    ncols = get_lpi_col_num(lpi)
    @assert(length(alpha) == ncols)

    
    if !isempty(row_name) 
        temp = Vector{String}()
        push!(temp,row_name)
        row_name = temp
    else 
        row_name = C_NULL 
    end

    #row_name = C_NULL

    beg = Ref{Cint}(0)
    
    ind = Vector{Cint}()
    entries = Vector{SCIP.SCIP_Real}()

    for i=1:ncols
        if abs(alpha[i]) < EPSILON continue end
        push!(ind,i-1)
        push!(entries,SCIP.SCIP_Real(alpha[i]))
    end
    
    if length(ind) >= 1
        SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(lpi, 1, lhs, rhs, row_name, length(ind), beg, ind, entries)
    end

end

function solve_lpi(
    lpi::Ptr{SCIP.SCIP_LPi}
)
    # Solve The LP
    if SCIP.SCIPlpiHasPrimalSolve() == 0
        @warn "LP Solver Does Not Support Primal Solve"
        throw(ErrorException("LP Solver Does Not Support Primal Solve"))
    end

    SCIP.@SCIP_CALL SCIP.SCIPlpiSolvePrimal(lpi)
    
    if SCIP.SCIPlpiIsPrimalInfeasible(lpi) == 1
        @warn "LP Solver Failed TO Solve Cut Generation LP: LP Infeasible"
        throw(ErrorException("LP Solver Does Not Support Primal Solve"))
    end
    if SCIP.SCIPlpiIsPrimalUnbounded(lpi) == 1
        @warn "LP Solver Failed TO Solve Cut Generation LP: LP Unbounded"
        throw(ErrorException("LP Solver Does Not Support Primal Solve"))
    end
    if SCIP.SCIPlpiIsOptimal(lpi) == 0
        @warn "LP Solver Failed TO Solve Cut Generation LP: LP Not optimzal" 
        throw(ErrorException("LP Solver Does Not Support Primal Solve"))
    end
end


"""
Get the LP solution from the LPI
"""
function get_lpi_solution_vector(lpi::Ptr{SCIP.SCIP_LPI})
    @assert(SCIP.SCIPlpiIsOptimal(lpi) != 0)

    dim = get_lpi_col_num(lpi)
    obj = Ref{SCIP.SCIP_Real}(0.0)
    sol = zeros(SCIP.SCIP_Real, dim)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(lpi, obj, sol,C_NULL, C_NULL, C_NULL)

    return sol
end

"""
Get the number of ccolumn in LPI
"""
function get_lpi_col_num(lpi::Ptr{SCIP.SCIP_LPI})::Int64
    ncols = Ref{Cint}(0)
    SCIP.SCIPlpiGetNCols(lpi,ncols)
    ncols = Int64(ncols[])
    return ncols
end