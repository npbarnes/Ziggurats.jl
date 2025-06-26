using Test, Ziggurats
using Random
using StatsBase
using Distributions
using SpecialFunctions

include("testutils.jl")

# Distributions for testing edge cases.
include("SteppedExponential.jl")
include("Doorstop.jl")

@testset "Ziggurats" begin
    # Test the test-distributions
    @testset "Stepped Exponential Distributions" begin
        test_distr(SteppedExponential(), 10^6)
        test_distr(SteppedExponential(0.5), 10^6)
        test_distr(SteppedExponential(1.5), 10^6)
    end

    @testset "Doorstop Distribution" begin
        test_distr(Doorstop(-1, 1, 3), 10^6)
        test_distr(Doorstop(-1, -1, 3), 10^6)
        test_distr(Doorstop(-1, 1, 1), 10^6)
    end

    # Test Ziggurats
    include("test_inverses.jl")
    include("argument_handling_tests.jl")
    include("interface_tests.jl")
    include("ziggurat_tests.jl")
end

nothing
