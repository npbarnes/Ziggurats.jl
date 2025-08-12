using Test, Ziggurats
using Random
using StatsBase
using Distributions
using SpecialFunctions
using AliasTables
using Aqua
using JET

include("init_tests.jl")

include("test_testdistributions.jl")

@testset "Aqua.jl Tests" begin
    Aqua.test_all(Ziggurats)
end

include("JET_tests.jl")

# Ziggurat tests
include("test_inverses.jl")
include("argument_handling_tests.jl")
include("interface_tests.jl")
include("completed_ziggurat_tests.jl")
include("sampling_tests.jl")

nothing
