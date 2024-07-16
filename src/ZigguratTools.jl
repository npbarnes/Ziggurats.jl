module ZigguratTools

using Roots
using Random
using Distributions
using Plots

export
    # Sampler types
    UnboundedIncreasingZiggurat,
    UnboundedDecreasingZiggurat,

    # Inverse pdfs
    ipdf_left,
    ipdf_right,

    # Plot ziggurats
    plotziggurat, plotziggurat!

abstract type AbstractZiggurat{X} <: Sampleable{Univariate, Continuous} end
Base.eltype(::AbstractZiggurat{X}) where X = X

include("inverse_pdf.jl")
include("unboundedmonotonic.jl")
include("boundedmonotonic.jl")
include("plotting.jl")

end # module
