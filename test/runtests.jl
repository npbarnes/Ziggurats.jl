using Test, ZigguratTools
using Random
using StatsBase

# Extension
DistributionsExt = Base.get_extension(ZigguratTools, :DistributionsExt)
@test DistributionsExt === nothing

using Distributions

DistributionsExt = Base.get_extension(ZigguratTools, :DistributionsExt)
@test DistributionsExt !== nothing

include("testutils.jl")

# Distributions for testing edge cases.
include("SteppedExponential.jl")
include("Doorstop.jl")

@testset "ZigguratTools" begin
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

    # Test ZigguratTools
    include("test_inverses.jl")
    include("ziggurat_tests.jl")
end

nothing
