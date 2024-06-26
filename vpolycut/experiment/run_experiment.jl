#
# run_experiment.jl
# The code to be called for an experiment run
#

module Experiment

# Code for parsing arguments
include("parseargs.jl")

# Code for the actual execution
include("execute.jl")

function main()
    settings = parseargs()
    execute(settings)
end

end

Experiment.main()