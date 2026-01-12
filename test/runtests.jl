include("init_tests.jl")

@testset "Aqua.jl Tests" begin
    Aqua.test_all(Ziggurats)
end

include("JET_tests.jl")

# Ziggurat tests
include("test_inverses.jl")
include("layermask_tests.jl")

include("layer_bits_tests.jl")

include("argument_handling_tests.jl")
include("assumptions_tests.jl")
include("interface_tests.jl")
include("test_testdistributions.jl")
include("completed_ziggurat_tests.jl")
include("sampling_tests.jl")

nothing
