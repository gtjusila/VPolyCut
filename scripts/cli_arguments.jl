using SCIP
using VPolyhedralCut
using ArgParse

struct ExecutionParameters
    instance::String
    separator_label::String
    separator::Type{<:SCIP.AbstractSeparator}
    easy::Bool
end

function get_execution_parameters()::ExecutionParameters
    cli_arguments = setup_cli_arguments()
    return read_cli_arguments(cli_arguments)
end

function setup_cli_arguments()
    settings = ArgParseSettings()
    @add_arg_table settings begin
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
    return parse_args(settings)
end

function read_cli_arguments(commandline_arguments::Dict)::ExecutionParameters
    instance = remove_whitespaces(commandline_arguments["instance"])
    mode_text = remove_whitespaces(commandline_arguments["mode"])
    separator = get_separator_type_from_string(mode_text)
    easy = commandline_arguments["easy"]
    return ExecutionParameters(instance, mode_text, separator, easy)
end

function get_separator_type_from_string(
    separator_type::String
)::Type{<:SCIP.AbstractSeparator}
    separator_type = remove_whitespaces(separator_type)
    if separator_type == "gomory"
        return GomorySeparator
    elseif separator_type == "vpc"
        return VPolyhedralCut.VPolyhedralSeparator
    else
        @error "Argument `mode` does not match any known modes. Using gomory mode"
        return GomorySeparator
    end
end

function remove_whitespaces(original_string::String)::String
    return string(strip(original_string))
end
