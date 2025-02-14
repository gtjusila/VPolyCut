module VPolyhedralCut

# SCIPJLUtils
include("scipjlutils/SCIPJLUtils.jl")

# Main VPolyhedralCut Module

# Tableau Object
# Handles All low level interaction with SCIP
# Also contains some abstract data structure which act as wrappers around SCIP data structures
include("tableauhandlers/Variable.jl")
include("tableauhandlers/LPRow.jl")
include("tableauhandlers/LPColumn.jl")
include("tableauhandlers/ConstraintMatrix.jl")
include("tableauhandlers/Tableau.jl")
include("tableauhandlers/scip_connector.jl")

# Common Data Structures
include("datastructures/Exceptions.jl")
include("datastructures/Point.jl")
include("datastructures/Ray.jl")
include("datastructures/NonBasicSpace.jl")
include("datastructures/CutPool.jl")
include("datastructures/PointRayCollection.jl")
include("datastructures/CornerPolyhedron.jl")

# Utilities
include("utilities/CPointer.jl")
include("utilities/scip_tableau_utilities.jl")
include("utilities/sepa_row_helpers.jl")
include("utilities/log_helpers.jl")
include("utilities/numerical_methods.jl")
include("utilities/eliminate_duplicate_rows.jl")
include("utilities/indicator_connector.jl")

# PRLP
include("prlp/PRLPstructures.jl")

# Branch and Bound
include("branchandbound/Node.jl")
include("branchandbound/NodeQueue.jl")
include("branchandbound/BranchingRule.jl")
include("branchandbound/BranchAndBound.jl")
include("branchandbound/scip_connector.jl")
include("branchandbound/execute.jl")

# Disjunction
include("disjunction/structures.jl")
include("disjunction/branch_and_bound_adaptor.jl")

# Datastructures for VPolyhedralCut
include("VPCStructures.jl")

# Parts 
include("subroutines/collect_point_rays.jl")
include("subroutines/construct_prlp.jl")
include("subroutines/gather_separating_solutions.jl")
include("subroutines/get_cut_from_separating_solution.jl")

# Separators
include("IndicatorSeparator.jl")
include("IntersectionSeparator.jl")
include("VPCexeclp.jl")

end