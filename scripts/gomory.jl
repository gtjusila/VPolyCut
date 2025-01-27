using ArgParse
using SCIP
using JuMP
using JSON
using VPolyhedralCut.SCIPJLUtils

function main()
    ### Command Line Argument ###
    args_setting = ArgParseSettings(;
        description="Solve the LP relaxation of the instance, do 1 round of gomory cuts and get the gap closed information. Output will be a results.json file containing gap information and scip solving statistic"
    )
    @add_arg_table args_setting begin
        "--instance", "-i"
        help = "An instance to be solve"
        required = true
        "--output_dir", "-o"
        help = "A directory where to writ results.json and scip_solving_statistic.txt"
        required = true
    end
    parameter = ArgParse.parse_args(args_setting)
    println(parameter)

    # Use VPolyhedralCut.SCIPJLUtils to create SCIP Optimizer object
    model = setup_scip_safe_jump_model()
    scip = get_scip_data_from_model(model)

    # Setup SCIP object parameter
    set_heuristics_emphasis_off(model)
    set_separators_emphasis_off(model)
    set_cut_selection_off(model)
    set_strong_branching_lookahead_off(model)

    JuMP.set_attribute(model, "limits/restarts", 0)
    JuMP.set_attribute(model, "limits/nodes", 1)
    JuMP.set_attribute(model, "limits/time", 3600) # Time Limit is usually not an issue for gomory runs
    JuMP.set_attribute(model, "separating/maxroundsroot", 1)
    #JuMP.set_attribute(model, "display/verblevel", 0)

    # Turn on gomory cut
    JuMP.set_attribute(model, "separating/gmi/freq", 0)
    JuMP.set_attribute(model, "separating/gmi/priority", 9999)
    JuMP.set_attribute(model, "separating/gmi/maxsuppabs", 5000)
    JuMP.set_attribute(model, "separating/gmi/dynamiccuts", false)
    JuMP.set_attribute(model, "separating/gmi/maxsupprel", 1.0)
    JuMP.set_attribute(model, "separating/gmi/forcecuts", true)

    # Read Problem
    instance_path = abspath(parameter["instance"])
    @info "Loading instance from $instance_path"
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(scip, instance_path, C_NULL)

    # Solve
    @info "Starting solve"
    SCIP.@SCIP_CALL SCIP.SCIPsolve(scip)

    # Write Output
    @info "Writing Result"
    result = Dict{String,Any}()
    result["separator"] = "gomory"
    result["instance"] = instance_path
    result["scip_status"] = SCIP.SCIPgetStatus(scip)
    result["initial_lp_obj"] = SCIP.SCIPgetFirstLPDualboundRoot(scip)
    result["final_lp_obj"] = SCIP.SCIPgetDualboundRoot(scip)

    output_path = abspath(parameter["output_dir"])
    result_path = joinpath(output_path, "results.json")
    open(result_path, "w") do io
        JSON.print(io, result, 4)
    end

    solving_statistic_path = joinpath(output_path, "scip_solving_statistic.txt")
    open(solving_statistic_path, "w") do io
        solv_stats_file_ptr = Libc.FILE(io)
        SCIP.SCIPprintStatistics(scip, solv_stats_file_ptr)
    end
end

main()