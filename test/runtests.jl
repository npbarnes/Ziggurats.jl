using SafeTestsets

@safetestset "Ziggurat Tests" begin
    @safetestset "Completed Ziggurat" include("completed_ziggurat_tests.jl")
    @safetestset "Sampling" include("sampling_tests.jl")
    @safetestset "Regression" include("regression_tests.jl")
end