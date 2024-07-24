module VPolyCut
include("constants.jl")

include("tableau/Variable.jl")
include("tableau/LPRow.jl")
include("tableau/LPColumn.jl")
include("tableau/ConstraintMatrix.jl")
include("tableau/Tableau.jl")
include("tableau/scip_connector.jl")

include("utils.jl")
include("CornerPolyhedron.jl")
include("IntersectionSeparator.jl")
end # module VPolyCut
