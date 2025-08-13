using Test, Ziggurats
using Random
using StatsBase
using Distributions
using SpecialFunctions
using AliasTables
using Aqua
using JET

include("stat_testutils.jl")

# Distributions for testing edge cases.
include("SteppedExponential.jl")
include("Doorstop.jl")

include("testdata.jl")

nothing
