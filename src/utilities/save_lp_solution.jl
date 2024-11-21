"""
    save_lp_solution(scip, lp_solution_path)

Save the current LP solution in `scip` to the file `lp_solution_path`
"""
function save_lp_solution(scip::SCIP.SCIPData, lp_solution_path::String)
    # Get the LP solution
    lp_solution = Ref{Ptr{SCIP.SCIP_Sol}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateLPSol(scip, lp_solution, C_NULL)

    # Write the LP solution to file
    solution_stream = open(lp_solution_path, "w")
    solution_stream_pointer = Libc.FILE(solution_stream)
    SCIP.@SCIP_CALL SCIP.SCIPprintSol(scip, lp_solution[], solution_stream_pointer, 0)

    # Free the memory
    SCIP.@SCIP_CALL SCIP.SCIPfreeSol(scip, lp_solution)

    @debug "LP Solution written to $(lp_solution_path)"
end