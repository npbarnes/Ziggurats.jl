module ZigguratTools

using Roots
using Random
using QuadGK
using LinearAlgebra
using AliasTables
using IntervalRootFinding
using ForwardDiff
using FixedPointNumbers
using Logging

export Ziggurat, MonotonicZiggurat
export BoundedZiggurat, UnboundedZiggurat, SymmetricZiggurat, CompositeZiggurat
export ziggurat, monotonic_ziggurat
export monotonic_segments, inversepdf

include("utilities.jl")
include("inverse_pdf.jl")
include("monotonic.jl")
include("composite.jl")

end # module
