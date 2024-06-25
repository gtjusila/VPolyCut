#
# writeoutput.jl
# Code to write output
#
import SCIP
function writeoutput(path::String, scip::SCIP.SCIPData, settings::Dict, reference_obj::SCIP.SCIP_Real)
    output = open(joinpath(path, "results.txt"), "w")

    println(output, "ExperimentMode: $(settings["mode"])")
    println(output, "Instance: $(settings["instance"])")
    println(output, "OutputDirectory: $(path)")

    println(output, "ReferenceObjective: $(round(reference_obj,sigdigits = 6))")

    initialdualbound = SCIP.SCIPgetFirstLPDualboundRoot(scip)
    println(output, "FirstLPObj: $(round(initialdualbound,sigdigits = 6))")

    initialgap = abs(initialdualbound - reference_obj)
    println(output, "InitialGap: $(round(initialgap,sigdigits = 6))")

    finaldualbound = SCIP.SCIPgetDualboundRoot(scip)
    println(output, "FinalLPObj: $(round(finaldualbound,sigdigits = 6))")

    finalgap = abs(finaldualbound - reference_obj)
    println(output, "FinalGap: $(round(finalgap, sigdigits=6 ))")

    gapclosed = 1 - finalgap / initialgap
    println(output, "GapClosed: $(round(gapclosed, sigdigits=6))")

    println(output, "SCIPStatus: $(SCIP.SCIPgetStatus(scip))")
    solv_stats = open(joinpath(path, "solving_statistic.txt"), "w")
    solv_stats_ptr = Libc.FILE(solv_stats)
    SCIP.SCIPprintStatistics(scip, solv_stats_ptr)

end