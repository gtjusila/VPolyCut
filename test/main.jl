# Setup Temp Forlders For Experiment
rm(joinpath("temp"), force=true, recursive=true)
mkdir(joinpath(pwd(),"temp"))

#include("./simple/simplex.jl")
#include("./soplex/soplex.jl")
include("./miplib/miplib.jl")