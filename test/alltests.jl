@testset "Ziggurat Tests" begin
    @testset "Completed Ziggurat" include("completed_ziggurat_tests.jl")
    @testset "Sampling" include("sampling_tests.jl")
end

nothing
