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
include("layermask_tests.jl")

# overlapping_bits_tests uses Supposition.jl which does not support x86.
# Future versions of Supposition may support x86. See Supposition.jl's issue #76 on GitHub
@static if Sys.ARCH !== :x86
    include("layer_bits_tests.jl")
end

include("argument_handling_tests.jl")
include("interface_tests.jl")
include("completed_ziggurat_tests.jl")
include("sampling_tests.jl")

nothing
