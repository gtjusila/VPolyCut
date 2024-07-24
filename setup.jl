using Pkg
using TOML

# Load configuration from TOML file
config = TOML.parsefile("config.toml")

# Extract SCIP.jl and SCIP directories from the configuration
scip_jl = strip(config["scip_jl"])
scip_dir = strip(config["scip_dir"])

# Remove existing Manifest.toml files
rm("./VPolyCut/Manifest.toml", force=true)
rm("./ExperimentScript/Manifest.toml", force=true)

# Setup for VPolyCut
cd("./vpolycut")
Pkg.activate(".")
ENV["SCIPOPTDIR"] = scip_dir
Pkg.develop(path=scip_jl)
Pkg.build()
try
    Pkg.test()
catch e
    println("Tests failed: ", e)
    exit(1) # Exit with an error code
end
Pkg.activate()

# Return to the root directory
cd("..")

# Setup for ExperimentScript
cd("./ExperimentScript")
Pkg.activate(".")
ENV["SCIPOPTDIR"] = scip_dir
Pkg.develop(path="../vpolycut")
Pkg.develop(path=scip_jl)
Pkg.build()
try
    Pkg.test()
catch e
    println("Tests failed: ", e)
    exit(1) # Exit with an error code
end
Pkg.activate()
