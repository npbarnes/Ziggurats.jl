using SafeTestsets

@safetestset "Ziggurat Tests" begin
    @safetestset "Completed Ziggurat" include("completed_ziggurat_tests.jl")
    @safetestset "Sampling Ziggurat" include("sampling_tests.jl")
end