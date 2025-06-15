module ZigguratTools

using Roots
using Random
using QuadGK
using LinearAlgebra
using AliasTables
using IntervalRootFinding
using ForwardDiff
using FixedPointNumbers

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

include("utilities.jl")
include("inverse_pdf.jl")
include("monotonic.jl")
include("composite.jl")

end # module
