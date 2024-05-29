import SCIP
include("constants.jl")

"""
Take a vector alpha of the problem dimension and add row lhs <= alpha*x <= rhs to LPI
"""
function add_row_to_lp(
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
    
    ncols = Ref{Cint}(0)
    SCIP.SCIPlpiGetNCols(lpi,ncols)
    ncols = Int64(ncols[])
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