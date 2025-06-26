module Ziggurats

using Roots
using Random
using QuadGK
using LinearAlgebra
using AliasTables
using IntervalRootFinding
using ForwardDiff
using FixedPointNumbers

export Ziggurat, MonotonicZiggurat
export BoundedZiggurat, UnboundedZiggurat, CompositeZiggurat
export ziggurat, monotonic_ziggurat
export inversepdf

public monotonic_segments, NoWrap

abstract type Ziggurat{X} end

include("utilities.jl")
include("inverse_pdf.jl")
include("monotonic.jl")
include("composite.jl")

function ziggurat(
    pdf,
    domain,
    N = nothing;
    ipdf = nothing,
    tailarea = nothing,
    cdf = nothing,
    ccdf = nothing,
    left_fallback = nothing,
    right_fallback = nothing,
    fallback = nothing,
    p = nothing
)
    domain = regularize(domain)

    if N === nothing
        N = eltype(domain) == Float16 ? 64 : 256
    end

    # regularize guarentees that length(domain) >= 2
    if length(domain) == 2
        monotonic_ziggurat(pdf, domain, N; ipdf, tailarea, cdf, ccdf, fallback)
    else
        CompositeZiggurat(
            pdf,
            domain,
            N;
            ipdfs = ipdf,
            cdf,
            ccdf,
            left_fallback,
            right_fallback,
            p
        )
    end
end

end # module
