struct CompositeZiggurat{X,Y,Z<:Tuple{Ziggurat{X,Y},Ziggurat{X,Y},Vararg{Ziggurat{X,Y}}},AT} <: Ziggurat{X,Y}
    zigs::Z
    at::AT
    function CompositeZiggurat(
        zigs::Tuple{Ziggurat{X,Y},Ziggurat{X,Y},Vararg{Ziggurat{X,Y}}},
        at::AliasTable
    ) where {X,Y}
        new{X,Y,typeof(zigs),typeof(at)}(zigs, at)
    end
end

# TODO: add option to autodetect monotonic subdomains.
function CompositeZiggurat(pdf, domain; kwargs...)
    domain = regularize(domain)
    N = default_numlayers(nothing, eltype(domain))
    CompositeZiggurat(pdf, domain, N; kwargs...)
end

function CompositeZiggurat(pdf, domain, N::Integer; kwargs...)
    domain = regularize(domain)
    Ns = fill(N, length(domain)-1)
    CompositeZiggurat(pdf, domain, Ns; kwargs...)
end

"""
    CompositeZiggurat(pdf, domain, Ns; [ipdf, cdf, ccdf, left_fallback, right_fallback, p])

Constructs a sampler for a piecewise-monotonic univariate probability distribution defined by
a probability density function (`pdf`). The domain must be a list of numbers that divides
the `pdf` into monotonic segments. The pdf must be monotonic on each subinterval. It must
not diverge to infinity anywhere on the domain, including at the division points, but may
otherwise be arbitrary - including discontinuous functions. Generate random numbers by
passing the returned ziggurat object to Julia's `rand` or `rand!` functions.

`CompositeZiggurat` independently generates monotonic ziggurats for each subdomain. The
leftmost and rightmost ziggurats may be unbounded. Internal subdomains are necessarily bounded.
The [`monotonic_ziggurat`](@ref) function is used internally. The `ipdf` argument is used to
pass a list of inverse pdfs; one for each subdomain. Use the `cdf` and `ccdf` arguments to
override the left and right tail areas respectively. Similarly, use `left_fallback` and
`right_fallback` if necessary.

# Arguments
 - `pdf`: The probability density function of the desired distribution. It must not diverge \
 to infinity anywhere on the `domain`, including the endpoints. It does not need to be \
 normalized, but `ipdf`, `cdf`, and `ccdf` need to have the same normalization as `pdf`.

 - `domain`: a list of numbers that divides the `pdf` into monotonic segments. `domain` may \
 The values will be promoted to the highest float type present, or Float64 if there \
 are no floats. The type that they are promoted to will be the type produced by sampling \
 from the resulting ziggurat. The domain may be unbounded in either or both directions, \
 for example, `(-Inf, 0, Inf)`.

 - `N`: (Optional) The number of layers in each ziggurat. `N` may be a single number that \
 applies to all ziggurats, or a list of numbers one for each subdomain. If `N` is a power \
 of two and the domain is Float64 with N <= 4096, Float32 with N <= 512, or Float16 with \
 N <= 64, then an optimized sampling algorithm is used. Normally, `N` defaults to 256, \
 but for Float16 domains it defaults to 64.

 - `ipdf`: (Optional) A list of functions that invert the pdf on each monotonic subdomain. \
 They are used in the ziggurat construction algorithm, but not during sampling, so they \
 may affect the performance of the initial call to `z = CompositeZiggurat(...)`, but not \
 subsequent calls to `rand(z)`. If `ipdf` is `nothing` or not provided, then all of the \
 ipdfs will be calculated numerically using a root finding algorithm. If any of the entries \
 in the list are `nothing` then those ipdfs will be calculated numerically. See also \
 [`inversepdf`](@ref).

 - `cdf`: (Optional) The cumulative distribution function. It should be normalized the same \
 way the pdf is normalized. Only used if the leftmost subdomain is unbounded below.

 - `ccdf`: (Optional) The complementary cumulative distribution function. It should be \
 normalized the same way the pdf is normalized. Only used if the rightmost subdomain is \
 unbounded above.

 - `left_fallback` and `right_fallback: (Optional) functions that take arguments (rng, a) \
 and produce a random variate in the tail beyond `a` (below `a` for `left_fallback`, and \
 above `a` for `right_fallback`). These are only used on unbounded domains, one or both \
 may be ignored if the domain is bounded on one or both sides. If no fallback is provided, \
 then inverse transform sampling will be used by numerically inverting the `cdf` and/or \
 the `ccdf` (which may themselves be numerical estimates). Overriding the fallbacks can \
 improve the worst-case performance of the sampling algorithm, but since the fallback is \
 only called rarely, it will have a smaller influence on the average performance.

# Examples
```julia-repl
julia> z = CompositeZiggurat(x->exp(-x^2), (-Inf, 0, Inf), 512)
CompositeZiggurat{...}(...)

julia> rand(z, 3)
3-element Vector{Float64}:
  0.756124133407901
 -1.014527969936336
 -1.3362254231230777
```
"""
function CompositeZiggurat(
    pdf,
    domain,
    Ns;
    ipdfs = nothing,
    cdf = nothing,
    ccdf = nothing,
    left_fallback = nothing,
    right_fallback = nothing,
    p = nothing
)
    domain = regularize(domain)

    a, b = extrema(domain)
    subdomains = get_subdomains(pdf, domain)

    if size(subdomains, 1) != length(Ns)
        throw(ArgumentError("Ns must be either an Integer or an iterable with length equal to the number of subdomains."))
    end
    if ipdfs === nothing
        ipdfs = [nothing for d in eachrow(subdomains)]
    end
    if cdf === nothing
        val = pdf(Roots.__middle(a, b))
        domain_type = eltype(domain)
        range_type = typeof(val)
        error_type = typeof(norm(val))
        segbuf = alloc_segbuf(domain_type, range_type, error_type)

        cdf = TailArea(pdf, a, segbuf)
    end
    if ccdf === nothing
        val = pdf(Roots.__middle(a, b))
        domain_type = typeof(domain[1])
        range_type = typeof(val)
        error_type = typeof(norm(val))
        segbuf = alloc_segbuf(domain_type, range_type, error_type)

        ccdf = TailArea(pdf, b, segbuf)
    end

    _p = zigprobs(cdf, subdomains, p)

    if !(length(domain) - 1 == length(Ns) == length(ipdfs) == length(_p))
        throw(ArgumentError("Ns, ipdfs, and p must all have a length one less than the length of domain."))
    end

    zig_gen = i -> begin
        if i == 1
            fallback = left_fallback
        elseif i == size(subdomains, 1)
            fallback = right_fallback
        else
            fallback = nothing
        end
        monotonic_ziggurat(pdf, @view(subdomains[i, :]), Ns[i]; ipdf = ipdfs[i], cdf, ccdf, fallback)
    end

    zigs = ntuple(zig_gen, size(subdomains, 1))
    at = AliasTable(_p)

    CompositeZiggurat(zigs, at)
end

"""
    monotonic_segments(f, domain)

Divides the domain into its monotonic subdomains using automatic differentiation and root
finding. The result will have every entry of `domain` as well as the detected critical points
of `f` in sorted order with duplicates dropped.
"""
function monotonic_segments(f, domain)
    domain = regularize(domain)
    df = x -> ForwardDiff.derivative(f, x)
    boundaries = find_zeros(df, domain)
    unique(sort!(vcat(domain.a, boundaries)))
end

"""
    get_subdomains(f, domain)

Convert a list of subdomain boundaries -- `domain` -- into a list of pairs representing subdomains.
Attempts to detect discontinuities in f, and account for them.

# Examples

```julia-repl
julia> Ziggurats.get_subdomains(cos, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, 0.0]
 [0.0, 1.0]

julia> Ziggurats.get_subdomains(x -> x>0 ? 1.0 : 0.0, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, 0.0]
 [5.0e-324, 1.0]

julia> Ziggurats.get_subdomains(x -> x>=0 ? 1.0 : 0.0, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, -5.0e-324]
 [0.0, 1.0]

julia> Ziggurats.get_subdomains(sign, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, -5.0e-324]
 [5.0e-324, 1.0]
```
"""
get_subdomains(f, domain) = _get_subdomains(f, regularize(domain))
get_subdomains(f, domain::Regularized) = _get_subdomains(f, domain)
function _get_subdomains(f, domain)
    sd = Array{eltype(domain)}(undef, length(domain)-1, 2)

    sd[1, 1] = domain[1]
    for (i, x) in enumerate(domain[2:(end - 1)])
        if isapprox(f(prevfloat(x)), f(x); atol = sqrt(eps(x)))
            sd[i, 2] = x
        else
            sd[i, 2] = prevfloat(x)
        end

        if isapprox(f(nextfloat(x)), f(x); atol = sqrt(eps(x)))
            sd[i + 1, 1] = x
        else
            sd[i + 1, 1] = nextfloat(x)
        end
    end
    sd[end, 2] = domain[end]
    sd
end

zigprobs(cdf, subdomains, ::Nothing) = [cdf(b) - cdf(a) for (a, b) in eachrow(subdomains)]
function zigprobs(cdf, subdomains, p::AbstractArray)
    if length(p) != size(subdomains, 1)
        error("the length of p is not equal to the number of subdomains.")
    end

    _p = similar(p, eltype(subdomains))
    for i in eachindex(p)
        if p[i] === nothing
            a, b = @view(subdomains[i, :])
            _p[i] = cdf(b) - cdf(a)
        else
            _p[i] = p[i]
        end
    end
    _p
end

@inline function Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:CompositeZiggurat})
    cz = s[]
    zigs = cz.zigs
    i = rand(rng, cz.at)
    @inbounds z = zigs[i]
    rand(rng, z)
end

@inline function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:CompositeZiggurat{X}}
) where {X<:FloatXX}
    cz = s[]
    if length(A) < 7 # TODO: Tune this number
        for i in eachindex(A)
            @inbounds A[i] = rand(rng, s)
        end
    else
        T = corresponding_uint(X)
        GC.@preserve A rand!(rng, Random.UnsafeView{T}(pointer(A), length(A)))

        for i in eachindex(A)
            @inbounds r = reinterpret(T, A[i])

            j = rand(rng, cz.at)
            @inbounds z = cz.zigs[j]

            w = widths(z)
            k = layerratios(z)
            y = heights(z)
            mb = highside(z)
            am = lowside(z)
            pdf = density(z)
            fb = fallback(z)
            LM = mask(z)

            @inbounds A[i] = zigsample(rng, w, k, y, mb, am, pdf, fb, LM, r)
        end
    end
    A
end
