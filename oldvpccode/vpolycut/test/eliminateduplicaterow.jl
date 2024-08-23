using VPolyCut
using SCIP
using SCIPExperimentUtils
using Test
@testset "eliminateduplicaterow" begin
    test = setup_scip_safe_jump_model()
    scip = get_scip_data_from_model(test)
    matrix::Matrix{Float64} = [1 2 3; 2 4 6; 1 0 2]
    @test VPolyCut.get_unique_row_indices(matrix, scip) == [1, 3]

    rows::Vector{Vector{Float64}} = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 0.0, 2.0]]
    @test VPolyCut.get_unique_row_indices(rows, scip) == [1, 3]
end