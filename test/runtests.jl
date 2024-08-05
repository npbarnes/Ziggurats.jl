using Test, ZigguratTools

# Setup
using Distributions
using Random
using StatsBase

include("testutils.jl")

# Distributions for testing edge cases.
include("SteppedExponential.jl")
include("Doorstop.jl")

# Test the test-distributions
@testset "Stepped Exponential Distributions" begin
    test_distr(SteppedExponential(), 10^6)
    test_distr(SteppedExponential(0.5), 10^6)
end

@testset "Doorstop Distribution" begin
    test_distr(Doorstop(-1, 1, 3), 10^6)
end

# Test ZigguratTools
include("ziggurat_tests.jl")
