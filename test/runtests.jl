using SCIP
using Revise
using JuMP
using VPolyhedralCut
using VPolyhedralCut.SCIPJLUtils
using Test

@testset "tableau" begin
    include("./tableautest.jl")
end

@testset "scipjl_utils" begin
    include("./scipjl_utils.jl")
end

@testset "intersectioncut" begin
    include("./intersectioncut.jl")
end

@testset "branchandbound" begin
    include("./branchandbound.jl")
end

@testset "leaves" begin
    include("./leaves.jl")
end

@testset "eliminateduplicaterow" begin
    include("./eliminateduplicaterow.jl")
end

@testset "vpolyhedralcut" begin
    include("./vpolyhedralcut.jl")
end