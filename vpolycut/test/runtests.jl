using SCIP
using Revise
using JuMP
using VPolyCut
using Test
using SCIPExperimentUtils

@testset "tableau" begin
    include("./tableautest.jl")
end

@testset "intersectioncut" begin
    include("./intersectioncut.jl")
end