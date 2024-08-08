module VPolyCut
include("utilities.jl")
include("numerical_methods.jl")

# Tableau
include("tableau/Variable.jl")
include("tableau/LPRow.jl")
include("tableau/LPColumn.jl")
include("tableau/ConstraintMatrix.jl")
include("tableau/Tableau.jl")
include("tableau/scip_connector.jl")
include("tableau/ComplementedTableau.jl")

# Branch and Bound
include("branchandbound/Node.jl")
include("branchandbound/NodeQueue.jl")
include("branchandbound/BranchingRule.jl")
include("branchandbound/BranchAndBound.jl")
include("branchandbound/scip_connector.jl")
include("branchandbound/execute.jl")

include("Projection.jl")
include("CornerPolyhedron.jl")
include("IntersectionSeparator.jl")
include("VPolyhedralSeparator.jl")
end # module VPolyCut
