module ZigguratTools

using Roots
using Random
using Distributions
using Plots

export
    # Sampler types
    UnboundedIncreasingZiggurat,
    UnboundedDecreasingZiggurat,


    # Plot ziggurats
    plotziggurat, plotziggurat!

abstract type AbstractZiggurat{X} <: Sampleable{Univariate, Continuous} end
Base.eltype(::AbstractZiggurat{X}) where X = X

include("unboundedmonotonic.jl")
include("boundedmonotonic.jl")
include("plotting.jl")

end # module
