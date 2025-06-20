module ZigguratTools

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
export monotonic_ziggurat
export inversepdf

public monotonic_segments, NoWrap

abstract type Ziggurat{X} end

include("utilities.jl")
include("inverse_pdf.jl")
include("monotonic.jl")
include("composite.jl")

end # module
