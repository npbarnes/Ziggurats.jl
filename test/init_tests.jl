using Test, Ziggurats
using Random
using StatsBase
using Distributions
using SpecialFunctions
using LambertW
using AliasTables
using Aqua
using JET
using Supposition

include("stat_testutils.jl")

# Distributions for testing edge cases.
include("SteppedExponential.jl")
include("Doorstop.jl")

include("testdata.jl")

nothing
