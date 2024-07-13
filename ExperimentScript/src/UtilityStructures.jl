import SCIP
import VPolyCut

struct ExecutionParameters
    instance::String
    separator_label::String
    separator::Type{<:SCIP.AbstractSeparator}
    easy::Bool
end

mutable struct ExperimentStore
    scip::SCIP.SCIPData
    result_path::String
    separator::String
    instance::String
    reference_objective::SCIP.SCIP_Real
    feasible::Bool
end

function ExperimentStore(scip, result_path, separator, instance)
    return ExperimentStore(scip, result_path, separator, instance, 0, false)
end