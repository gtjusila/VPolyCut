struct TimeLimitExceeded <: Exception end
struct LPSolNotBasic <: Exception end
struct LPSolNotOptimal <: Exception end
struct FailedDisjunctiveLowerBoundTest <: Exception end
struct FailedToProvePRLPFeasibility <: Exception end