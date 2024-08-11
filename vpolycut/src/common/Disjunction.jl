# A special data structure for managing disjunctions
# It pair a point and row
@kwdef mutable struct Disjunction
    points::Vector{CornerPoint}
    row::Vector{SCIP_Real}
end