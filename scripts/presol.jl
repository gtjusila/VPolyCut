using ArgParse
using SCIP
using JuMP
using JSON
using TOML
using VPolyhedralCut.SCIPJLUtils

mutable struct FirstLPEvent <: SCIP.AbstractEventhdlr
    scip::SCIP.SCIPData
    n_frac::Int
    nlp_rows::Int
    nlp_cols::Int
end

function SCIP.eventexec(event::FirstLPEvent)
    @info "Starting Event"
    event.n_frac = SCIP.SCIPgetNLPBranchCands(event.scip)
    @info "Number of Fractional Branching Candidates: $(event.n_frac)"
    event.nlp_rows = SCIP.SCIPgetNLPRows(event.scip)
    @info "Number of LP Rows: $(event.nlp_rows)"
    event.nlp_cols = SCIP.SCIPgetNLPCols(event.scip)
    @info "Number of LP Columns: $(event.nlp_cols)"
end

function main()
    ### Command Line Argument ###
    args_setting = ArgParseSettings(;
        description = "Solve the LP relaxation of the instance, do 1 round of gomory cuts and get the gap closed information. Output will be a results.json file containing gap information and scip solving statistic"
    )
    @add_arg_table args_setting begin
        "--instance", "-i"
        help = "An instance to be solve"
        required = true
        "--output_dir", "-o"
        help = "A directory where to writ results.json and scip_solving_statistic.txt"
        required = true
        "--config", "-c"
        help = "Parameters for the VPolyhedralCut algorithm"
        required = true
        "--solution", "-s"
        help = "Solution"
        default = ""
    end
    parameter = ArgParse.parse_args(args_setting)

    # Setup output directory
    output_path = abspath(parameter["output_dir"])
    log_path = joinpath(output_path, "vpc_logs.log")
    config = TOML.parsefile(parameter["config"])

    # Use VPolyhedralCut.SCIPJLUtils to create SCIP Optimizer object
    model = setup_scip_safe_jump_model()
    scip = get_scip_data_from_model(model)

    # Setup SCIP object parameter
    if (!config["scip_enable_heuristic"])
        set_heuristics_emphasis_off(model)
    end
    if (!config["scip_enable_conflict_analysis"])
        JuMP.set_attribute(model, "conflict/enable", false)
    end
    if (config["scip_disable_scip_cuts"])
        set_separators_emphasis_off(model)
    end
    if (!config["scip_enable_cut_selection"])
        set_cut_selection_off(model)
    end
    if (!config["scip_enable_strong_branching_lookahead"])
        set_strong_branching_lookahead_off(model)
    end
    if (!config["scip_enable_root_node_propagation"])
        set_root_node_propagation_off(model)
    end
    if !config["scip_allow_restart"]
        JuMP.set_attribute(model, "limits/restarts", 0)
        JuMP.set_attribute(model, "estimation/restarts/restartpolicy", 'n')
        JuMP.set_attribute(model, "presolving/maxrestarts", 0)
    end
    JuMP.set_attribute(model, "limits/nodes",
        config["scip_node_limit"]
    )
    JuMP.set_attribute(model, "limits/time",
        config["scip_time_limit"]
    )
    JuMP.set_attribute(
        model, "separating/maxroundsroot", config["scip_max_root_cutting_plane_rounds"]
    )
    # Read Problem
    instance_path = abspath(parameter["instance"])
    @info "Loading instance from $instance_path"
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(scip, instance_path, C_NULL)
    eventhdlr = FirstLPEvent(scip, 0, 0, 0)
    SCIP.include_event_handler(
        scip,
        eventhdlr;
        name = "firstlp"
    )
    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(scip)
    SCIP.catch_event(scip, SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED, eventhdlr)
    # Solve
    @info "Starting solve"
    SCIP.@SCIP_CALL SCIP.SCIPsolve(scip)

    # Write Output
    @info "Writing Result"
    result = Dict{String,Any}()
    result["instance"] = instance_path
    result["n_frac"] = eventhdlr.n_frac
    result["nlp_rows"] = eventhdlr.nlp_rows
    result["nlp_cols"] = eventhdlr.nlp_cols

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