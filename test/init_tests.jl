using Test, Ziggurats
using Random
using StatsBase
using Distributions
using SpecialFunctions
using AliasTables
using Aqua
using JET

# Future versions of Supposition may support x86. See Supposition.jl's issue #76 on GitHub
@static if Sys.ARCH !== :x86
    using Supposition
end

include("stat_testutils.jl")

# Distributions for testing edge cases.
include("SteppedExponential.jl")
include("Doorstop.jl")

include("testdata.jl")

nothing
