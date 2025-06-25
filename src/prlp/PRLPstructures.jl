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
    BARRIER_WITH_CROSSOVER = 4
end
@enum TerminationStatus LPI_OPTIMAL LPI_TIME_LIMIT_EXCEEDED LPI_NOT_SOLVED LPI_UNBOUNDED LPI_INFEASIBLE

@kwdef mutable struct PRLP
    solve_algorithm::PRLPsolveAlgorithm = PRIMAL_SIMPLEX

    # Permanent Data
    scip::SCIP.SCIPData
    dimension::Int = 0
    points::Vector{Point} = []
    rays::Vector{Ray} = []
    lpi::CPtr{SCIP.SCIP_LPI} = CPtr(SCIP.SCIP_LPI)
    allow_warm_start::Bool = true
    beta::SCIP.SCIP_Real = 1.0

    # State
    lp_constructed::Bool = false
    last_solve_time::SCIP.SCIP_Real = 0.0
    last_simplex_iterations::Int = 0.0
    solution_available::Bool = false
    solution_vector::RealVector = []
    solution_objective::SCIP.SCIP_Real = 0.0
    last_solve_good::Bool = false
    last_good_state::Ref{Ptr{SCIP.SCIP_LPISTATE}} = Ref(C_NULL)

    # Statistic
    solve_statistics::Vector{Any} = []
    n_basis_restart::Int = 0
end