using SCIP
using Printf

"""
`RealVector` is an alias for `SCIP.SCIP_Real`
"""
const RealVector = Vector{SCIP.SCIP_Real}
# Base Struct for the point ray linear program

@enum PRLPsolveAlgorithm begin
    PRIMAL_SIMPLEX = 1
    DUAL_SIMPLEX = 2
    BARRIER = 3
end
@enum TerminationStatus LPI_OPTIMAL LPI_TIME_LIMIT_EXCEEDED LPI_NOT_SOLVED LPI_UNBOUNDED

@kwdef mutable struct PRLP
    dimension::Int = 0
    points::Vector{Point} = []
    rays::Vector{Ray} = []
    lp_constructed::Bool = false
    lpi::CPtr{SCIP.SCIP_LPI} = CPtr(SCIP.SCIP_LPI)
    last_solve_time::SCIP.SCIP_Real = 0.0
    last_simplex_iterations::Int = 0.0
    solution_available::Bool = false
    solution_vector::RealVector = []
    solution_objective::SCIP.SCIP_Real = 0.0
    solve_algorithm::PRLPsolveAlgorithm = PRIMAL_SIMPLEX
    solve_statistics::Vector{Any} = []
end