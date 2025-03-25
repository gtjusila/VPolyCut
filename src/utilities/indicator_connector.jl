using SCIP

"""
    IndicatorVarEntry

An object to represent an indicator variable in SCIP. This object contains the original binary variable,
whether there are constraints involving its negation, the constriants that is tied to the original binary variable, 
and the constraint that is tied to the negated binary variable along with their associated slacks.
"""
mutable struct IndicatorVarEntry
    indicator_var::Ptr{SCIP.SCIP_VAR} # We always use original var as key
    is_negated::Bool
    negated_var::Ptr{SCIP.SCIP_VAR}
    constraints::Vector{Ptr{SCIP.SCIP_CONS}} # Indicator constraints that are tied to the original binary variable
    negated_constraints::Vector{Ptr{SCIP.SCIP_CONS}} # Indicator constraints that are tied to the negated binary variable
end

function is_indicator_satisfied(scip::SCIP.SCIPData, entry::IndicatorVarEntry)::Bool
    for constraint in entry.constraints
        if is_true(SCIP.SCIPisViolatedIndicator(scip, constraint, C_NULL))
            return false
        end
    end
    for constraint in entry.negated_constraints
        if is_true(SCIP.SCIPisViolatedIndicator(scip, constraint, C_NULL))
            return false
        end
    end
    return true
end

function is_negated_variable(var::Ptr{SCIP.SCIP_VAR})::Bool
    return SCIP.SCIPvarGetStatus(var) == SCIP.SCIP_VARSTATUS_NEGATED
end

function getIndicatorConstraints(scip::SCIP.SCIPData)::Vector{Ptr{SCIP.SCIP_CONS}}
    indicator_constrainthandler = SCIP.SCIPfindConshdlr(scip, "indicator")
    n_indicator_constraints = SCIP.SCIPconshdlrGetNConss(indicator_constrainthandler)
    indicator_constraints_array = SCIP.SCIPconshdlrGetConss(indicator_constrainthandler)
    return unsafe_wrap(
        Vector{Ptr{SCIP.SCIP_CONS}}, indicator_constraints_array, n_indicator_constraints
    )
end

"""
    getIndicatorVarDict(scip::SCIP.SCIPData)

Given a SCIPData object. This function retrieves information about the indicator constraints in the SCIP Object.
It returns a dictionary with the original binary variable as the key and the IndicatorVarEntry as the value.
"""
function getIndicatorVarDict(
    scip::SCIP.SCIPData
)::Dict{Ptr{SCIP.SCIP_VAR},IndicatorVarEntry}
    indicator_constraints = getIndicatorConstraints(scip)
    indicator_var_dict = Dict{Ptr{SCIP.SCIP_VAR},IndicatorVarEntry}()

    for constraint in indicator_constraints
        orig_binary_variable = SCIP.SCIPgetBinaryVarIndicator(constraint)
        negated_binary_variable = C_NULL
        is_negated = false

        if is_negated_variable(orig_binary_variable)
            is_negated = true
            negated_binary_variable = orig_binary_variable
            orig_binary_variable = SCIP.SCIPvarGetNegationVar(orig_binary_variable)
        end

        if !haskey(indicator_var_dict, orig_binary_variable)
            indicator_var_dict[orig_binary_variable] = IndicatorVarEntry(
                orig_binary_variable, is_negated, negated_binary_variable, [], []
            )
        end

        if is_negated
            indicator_var_dict[orig_binary_variable].is_negated = true
            push!(
                indicator_var_dict[orig_binary_variable].negated_constraints,
                constraint
            )
        else
            push!(
                indicator_var_dict[orig_binary_variable].constraints,
                constraint
            )
        end
    end
    return indicator_var_dict
end
