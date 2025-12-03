module Ziggurats

using Roots
using Random
using QuadGK
using LinearAlgebra
using AliasTables
using IntervalRootFinding
using ForwardDiff
using FixedPointNumbers
using ILog2
using Logging
using Compat

export Ziggurat, MonotonicZiggurat
export BoundedZiggurat, UnboundedZiggurat, CompositeZiggurat
export ziggurat, monotonic_ziggurat
export inversepdf

@compat public monotonic_segments, NoWrap

abstract type Ziggurat{X} end

include("utilities.jl")
include("inverse_pdf.jl")
include("monotonic.jl")
include("composite.jl")

"""
    ziggurat(pdf, domain, [N]; [ipdf, tailarea, cdf, ccdf, fallback, left_fallback, right_fallback, p])

Constructs a high-performance sampler for a univariate probability distribution defined by a
probability density function (`pdf`). The domain must be a list of numbers that divides the
`pdf` into monotonic segments. The pdf must be monotonic on each subinterval. It must not
diverge to infinity anywhere on the domain, including at the division points, but may
otherwise be arbitrary - including discontinuous functions. Generate random numbers by
passing the returned ziggurat object to Julia's `rand` or `rand!` functions.

If the domain only has two endpoints and no internal points, then the keyword arguments
`ipdf`, `tailarea`, `cdf`, `ccdf`, and `fallback` get passed along to
[`monotonic_ziggurat`](@ref), otherwise the arguments  `ipdf`, `cdf`, `ccdf`, `left_fallback`,
`right_fallback`, and `p` get passed along to [`CompositeZiggurat`](@ref). See the
documentation of those functions for more details on the interactions between arguments.

# Examples
```julia-repl
julia> using Ziggurats
julia> z = ziggurat(x -> exp(-x^2/2), (-Inf, 0, Inf));
```
`Ziggurats.jl` hooks into Julia's `Random` API, so the various incarnations of `rand` and
`rand!` will work as expected.

```julia-repl
julia> rand(z)
0.6223352924397899

julia> rand(z, 3)
3-element Vector{Float64}:
 0.7997105858285334
 0.24692390755438537
 0.5731782397234146

julia> using Random

julia> a = Vector{Float64}(undef, 3);

julia> rand!(a, z)
3-element Vector{Float64}:
  1.861164263929902
  0.19287723622925435
 -0.21573227151484722
```
"""
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
        CompositeZiggurat(pdf, domain, N; ipdfs = ipdf, cdf, ccdf, left_fallback, right_fallback, p)
    end
end

end # module
