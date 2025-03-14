using ArgParse
using SCIP
using JuMP
using JSON
using TOML
using VPolyhedralCut
using VPolyhedralCut.SCIPJLUtils

function main()
    ### Command Line Argument ###
    args_setting = ArgParseSettings(;
        description = "Solve the LP relaxation of the instance, do 1 round of vpc cuts and get the gap closed information. Output will be a results.json file containing gap information and scip solving statistic"
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
    #JuMP.set_attribute(model, "display/verblevel", 0)

    # Turn on vpc cut
    vpcparam = VPolyhedralCut.VPCParameters(;
        n_leaves = config["vpc_n_leaves"],
        log_directory = output_path,
        time_limit = 900,
        min_restart = config["vpc_min_restart"],
        max_round = config["vpc_max_participating_round"],
        prlp_solve_method = config["vpc_prlp_solve_method"],
        prlp_allow_warm_start = config["vpc_prlp_allow_warm_start"],
        cut_limit = config["vpc_max_cut_per_round"],
        prlp_max_consecutive_fail = config["vpc_prlp_max_consecutive_fail"],
        prlp_min_increase_non_stagnating = config["vpc_prlp_min_gap_closed_increase"],
        prlp_max_consecutive_stagnation = config["vpc_prlp_max_consecutive_stagnation"]
    )

    vpcsepa = VPolyhedralCut.include_vpolyhedral_sepa(
        scip; parameters = vpcparam,
        delay = config["vpc_delayed"], priority = config["vpc_priority"],
        freq = config["vpc_frequency"])

    # Read Problem
    instance_path = abspath(parameter["instance"])
    @info "Loading instance from $instance_path"
    SCIP.@SCIP_CALL SCIP.SCIPreadProb(scip, instance_path, C_NULL)

    # If solution is given then read sol
    if parameter["solution"] != ""
        SCIP.@SCIP_CALL SCIP.SCIPreadSol(scip, parameter["solution"])
    end
    # Solve

    @info "Starting solve"
    SCIP.@SCIP_CALL SCIP.SCIPsolve(scip)

    # Write Output
    @info "Writing Result"
    result = Dict{String,Any}()
    result["separator"] = "vpc"
    result["instance"] = instance_path
    result["scip_status"] = SCIP.SCIPgetStatus(scip)
    result["initial_lp_obj"] = SCIP.SCIPgetFirstLPDualboundRoot(scip)
    result["final_lp_obj"] = SCIP.SCIPgetDualboundRoot(scip)
    result["node_count"] = SCIP.SCIPgetNNodes(scip)
    result["final_gap"] = SCIP.SCIPgetGap(scip)
    result["solve_time"] = SCIP.SCIPgetSolvingTime(scip)
    result["sepa_termination_message"] = vpcsepa.termination_status
    result["parameters"] = vpcparam
    result["statistics"] = vpcsepa.statistics
    result["num_runs"] = SCIP.SCIPgetNRuns(scip)
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