module ZigguratTools

using Roots
using Random
using QuadGK
using LinearAlgebra

export
    # Abstract types
    MonotonicZiggurat,

    # Sampler types
    UnboundedZiggurat,
    BoundedZiggurat

# Utility
function between(a, b, x)
    l, r = minmax(a, b)
    l <= x <= r
end

include("inverse_pdf.jl")
include("monotonic.jl")

end # module
