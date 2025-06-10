module ZigguratTools

using Roots
using Random
using QuadGK
using LinearAlgebra
using AliasTables
using IntervalRootFinding
using ForwardDiff

export
    # Abstract types
    MonotonicZiggurat,

    # Sampler types
    UnboundedZiggurat,
    BoundedZiggurat,
    SymmetricZiggurat,
    CompositeZiggurat,

    # Smart constructors
    monotonic_ziggurat,

    # Helper functions,
    monotonic_segments

public inverse

# Utility
function between(a, b, x)
    l, r = minmax(a, b)
    l <= x <= r
end

function regularize_domain(domain)
    unique(sort!(collect(promote(float.(domain)...))))
end

include("inverse_pdf.jl")
include("monotonic.jl")
include("composite.jl")

end # module
