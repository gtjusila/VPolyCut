module VPolyhedralCut

# SCIPJLUtils
include("scipjlutils/SCIPJLUtils.jl")

# Common Utilities
include("commons/CPointer.jl")
include("commons/numerical_methods.jl")
include("commons/scip_tableau_utilities.jl")
include("commons/log_helpers.jl")
include("commons/Exceptions.jl")

# Nonbasicspace (Point and Ray type is needed by Nonbasicspace)
include("pointray/Point.jl")
include("pointray/Ray.jl")
include("nonbasicspace/ConstraintMatrix.jl")
include("nonbasicspace/NonBasicSpace.jl")

# Branch and Bound
include("branchandbound/Node.jl")
include("branchandbound/NodeQueue.jl")
include("branchandbound/BranchingRule.jl")
include("branchandbound/BranchAndBound.jl")
include("branchandbound/scip_connector.jl")
include("branchandbound/execute.jl")

# Disjunction
include("disjunction/structures.jl")

# Point Ray Collection
include("pointray/CornerPolyhedron.jl")
include("pointray/PointRayCollection.jl")

# Utilities
include("utilities/indicator_connector.jl")

# PRLP
include("prlp/PRLPstructures.jl")
include("prlp/ObjectivePool.jl")
include("prlp/PRLPmethods.jl")

# CutPool
include("cutpool/CutPool.jl")
include("cutpool/sepa_row_helpers.jl")

# Datastructures for VPolyhedralCut
# and stuff that rely on VPC structures
include("VPCStructures.jl")
include("disjunction/branch_and_bound_adaptor.jl")
include("pointray/collect_point_rays.jl")
include("prlp/construct_prlp.jl")
include("prlp/gather_separating_solutions.jl")
include("cutpool/get_cut_from_separating_solution.jl")
include("utilities/estimate_prlp_nonzero.jl")
include("utilities/get_analytic_center.jl")
# Separators
#include("IndicatorSeparator.jl")
#include("IntersectionSeparator.jl")
include("VPCexeclp.jl")

end