using Pkg
using Logging: Logging
Pkg.activate(".")
using Coverage
try
    Pkg.test(; coverage=true)
catch e
    println(e)
end
Logging.disable_logging(Logging.Info)
coverage = process_folder()
LCOV.writefile("lcov.info", coverage)
covered, total = get_summary(coverage)
println("Covered lines: $covered. Total lines: $total. Coverage: $(covered/total*100)%")
Coverage.clean_folder(".")