module ZigguratTools

using Roots
using Random
using Distributions
using Plots

export
    # Abstract types
    Ziggurat,
    MonotonicZiggurat,

    # Sampler types
    UnboundedMonotonicZiggurat,
    BoundedMonotonicZiggurat,

    # Helper
    monotonic_ziggurat,

    # Inverse pdfs
    ipdf_left,
    ipdf_right,
    inverse,

    # Plot ziggurats
    plotziggurat,
    plotziggurat!

abstract type Ziggurat{X} <: Sampleable{Univariate,Continuous} end
Base.eltype(::Ziggurat{X}) where {X} = X

function between(a, b, x)
    l, r = minmax(a, b)
    l <= x <= r
end

include("inverse_pdf.jl")
include("monotonic.jl")
include("plotting.jl")

end # module
