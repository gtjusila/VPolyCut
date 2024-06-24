#
# parseargs.jl
# Contains the code for CLI parsing arguments
#
import ArgParse
"""
Parse Command Line Arguments
"""
function parseargs()
    settings = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table settings begin
        "--instance", "-i"
        help = "Name of instance. instance.sol and instance.mps should be in the data folder"
        required = true
        "--mode", "-m"
        help = "The mode of experiment to be run. Modes available: gomory, vpoly"
        required = true
    end
    return ArgParse.parse_args(settings)
end