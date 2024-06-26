#
# helpers.jl
# Various help funciton for running experiment
# 
import SCIP
import Dates
using Printf

function setupexperimentdirectory(instance::String, mode::String)
    cwd = pwd()
    results_dir = joinpath(cwd, "results")

    # Check if the results directory exists, if not, create it
    if !isdir(results_dir)
        mkpath(results_dir)
        println("Created directory: results")
    end

    # Get the current date and time
    now = Dates.now()
    timestamp = Dates.format(now, "yyyymmddTHHMMSS")

    # Create the new directory name
    new_dir_name = @sprintf("%s_%s_%s", timestamp, instance, mode)
    new_dir_path = joinpath(results_dir, new_dir_name)

    # Check if the directory already exists
    if isdir(new_dir_path)
        println("You are already running a similar experiment currently.")
        exit(1)
    else
        mkpath(new_dir_path)
        println("Created directory: $new_dir_path")
    end

    return new_dir_path
end

"""
Load solution from path and place it to pointer
"""
function load_solution(scip::SCIP.SCIPData, sol::Ref{Ptr{SCIP.SCIP_Sol}}, path::String)
    partial = Ref{SCIP.SCIP_Bool}(false)
    error = Ref{SCIP.SCIP_Bool}(false)
    SCIP.@SCIP_CALL SCIP.SCIPcreateOrigSol(scip, sol, C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPreadSolFile(scip, path, sol[], false, partial, error)
    if error[] == true
        @error "Error when reading solution"
    end
    if partial[] == true
        @error "Solution cannot a parial solution"
    end
end