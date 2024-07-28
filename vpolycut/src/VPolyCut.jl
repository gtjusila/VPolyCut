module VPolyCut
include("typedefinitions.jl")

# Tableau
include("tableau/Variable.jl")
include("tableau/LPRow.jl")
include("tableau/LPColumn.jl")
include("tableau/ConstraintMatrix.jl")
include("tableau/Tableau.jl")
include("tableau/scip_connector.jl")

# Branch and Bound
include("branchandbound/Node.jl")
include("branchandbound/BranchAndBound.jl")
include("branchandbound/scip_connector.jl")
include("branchandbound/execute.jl")

include("helper.jl")
include("Projection.jl")
include("CornerPolyhedron.jl")
include("IntersectionSeparator.jl")
end # module VPolyCut
