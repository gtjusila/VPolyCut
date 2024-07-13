import SCIP
import VPolyCut
import ArgParse

"""
get_execution_parameters

Get ExecutionParameters from commandline_arguments
"""
function get_execution_parameters()::ExecutionParameters
    commandline_arguments = get_commandline_arguments()
    return parse_commandline_arguments(commandline_arguments)
end

"""
get_commandline_arguments

Parse Command Line Arguments
"""
function get_commandline_arguments()
    settings = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table settings begin
        "--instance", "-i"
        help = "Name of instance. instance.sol and instance.mps should be in the data folder"
        required = true
        "--mode", "-m"
        help = "The mode of experiment to be run. Modes available: gomory, vpoly"
        required = true
        "--easy", "-e"
        help = "Disable Presolving and Propagation for easy instances"
        action = :store_true
    end
    return ArgParse.parse_args(settings)
end

"""
parse_commandline_arguments

Create an object of type ExecutionParameters from the given `commandline_arguments` dicitonary
"""
function parse_commandline_arguments(commandline_arguments::Dict)::ExecutionParameters
    instance = remove_whitespaces(commandline_arguments["instance"])
    mode_text = remove_whitespaces(commandline_arguments["mode"])
    separator = get_separator_type_from_string(mode_text)
    easy = commandline_arguments["easy"]
    return ExecutionParameters(instance, mode_text, separator, easy)
end

"""
get_seperator_type_from_string

take string `separator_type` and convert it to the coresponding seperator_type
"""
function get_separator_type_from_string(separator_type::String)::Type{<:SCIP.AbstractSeparator}
    separator_type = remove_whitespaces(separator_type)
    if separator_type == "gomory"
        return GomorySeparator
    elseif separator_type == "vpc"
        return VPolyCut.IntersectionSeparator
    else
        @error "Argument `mode` does not match any known modes. Using gomory mode"
        return GomorySeparator
    end
end

function remove_whitespaces(original_string::String)::String
    return string(strip(original_string))
end
