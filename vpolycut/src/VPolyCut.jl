module VPolyCut
include("numerical_methods.jl")

# Tableau Object
# Handles All low level interaction with SCIP
include("tableau/Variable.jl")
include("tableau/LPRow.jl")
include("tableau/LPColumn.jl")
include("tableau/ConstraintMatrix.jl")
include("tableau/Tableau.jl")
include("tableau/scip_connector.jl")
include("tableau/ComplementedTableau.jl")

# utilities
include("utilities.jl")

# Branch and Bound
include("branchandbound/Node.jl")
include("branchandbound/NodeQueue.jl")
include("branchandbound/BranchingRule.jl")
include("branchandbound/BranchAndBound.jl")
include("branchandbound/scip_connector.jl")
include("branchandbound/execute.jl")

# Common Data Structures
include("common/Projection.jl")
include("common/eliminate_duplicate.jl")
include("common/CutPool.jl")
include("common/PointRayCollection.jl")
include("common/CornerPolyhedron.jl")

include("IntersectionSeparator.jl")
include("VPolyhedralSeparator.jl")
end # module VPolyCut
