using VPolyhedralCut.SCIPJLUtils
using SCIP
using JuMP

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
    constraints::Vector{Pair{Ptr{SCIP.SCIP_CONS},Ptr{SCIP.SCIP_VAR}}}
    negated_constraints::Vector{Pair{Ptr{SCIP.SCIP_CONS},Ptr{SCIP.SCIP_VAR}}}
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
        slack_variable = SCIP.SCIPgetSlackVarIndicator(constraint)
        linear_constraint = SCIP.SCIPgetLinearConsIndicator(constraint)
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
                constraint => slack_variable
            )
        else
            push!(
                indicator_var_dict[orig_binary_variable].constraints,
                constraint => slack_variable
            )
        end
    end
    return indicator_var_dict
end

model = setup_scip_safe_jump_model()
set_everything_off(model)
scip = get_scip_data_from_model(model)

#JuMP.set_attribute(model, "limits/nodes", 1)
JuMP.set_attribute(model, "constraints/indicator/sepafreq", -1)

# Load problem
SCIP.@SCIP_CALL SCIP.SCIPreadProb(
    scip, "experiment_data/indicator_instances/supportcase21i.mps", C_NULL
)

my_function = function (data)
    indicator_var_dict = getIndicatorVarDict(data["scip"])
    for key in keys(indicator_var_dict)
        println("Var: ", unsafe_string(SCIP.SCIPvarGetName(key)))
        println("Is Negated: ", indicator_var_dict[key].is_negated)
        println(
            "Negated Var: ",
            unsafe_string(
                SCIP.SCIPvarGetName(indicator_var_dict[key].negated_var)
            )
        )
        println(
            "Constraints: ",
            [
                unsafe_string(SCIP.SCIPconsGetName(pair.first)) for
                pair in indicator_var_dict[key].constraints
            ]
        )
        println(
            "Negated Constraints: ",
            [
                unsafe_string(SCIP.SCIPconsGetName(pair.first)) for
                pair in indicator_var_dict[key].negated_constraints
            ]
        )
    end
end

data = Dict("scip" => scip)

my_callback = initiate_callback(
    my_function,
    model,
    data
)

SCIP.@SCIP_CALL SCIP.SCIPtransformProb(scip)

register_callback(
    model,
    SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED,
    my_callback
)

SCIP.@SCIP_CALL SCIP.SCIPsolve(scip)