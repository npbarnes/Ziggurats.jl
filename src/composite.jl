struct CompositeZiggurat{X,Y,Z<:BareZiggurat{X,Y},F,LF,RF,S,AT<:AliasTable} <: Ziggurat{X}
    zigs::Vector{Z}
    pdf::F
    left_fallback::LF
    right_fallback::RF
    shifts::Vector{S}
    at::AT
end

# TODO: add option to autodetect monotonic subdomains.
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
override the left and right tail areas respectively. Similarly, use `left_fallback` a and
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

    if length(subdomains) != length(Ns)
        throw(ArgumentError("Ns must be either an Integer or an iterable with length equal to the number of subdomains."))
    end
    if ipdfs === nothing
        ipdfs = [nothing for d in subdomains]
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
    if p === nothing
        _p = [cdf(b) - cdf(a) for (a, b) in subdomains]
    else
        _p = similar(p, Float64)
        for i in eachindex(p)
            if p[i] === nothing
                a, b = subdomains[i]
                _p[i] = cdf(b) - cdf(a)
            else
                _p[i] = p[i]
            end
        end
    end

    if !(length(domain) - 1 == length(Ns) == length(ipdfs) == length(_p))
        throw(ArgumentError("Ns, ipdfs, and p must all have a length one less than the length of domain."))
    end

    zig_gen = i -> begin
        if i == 1
            fallback = left_fallback
        elseif i == length(subdomains)
            fallback = right_fallback
        else
            fallback = nothing
        end
        monotonic_ziggurat(pdf, subdomains[i], Ns[i]; ipdf = ipdfs[i], cdf, ccdf, fallback)
    end

    zigs = [zig_gen(i) for i in 1:length(subdomains)]

    barezigs = bareziggurat.(zigs)
    shifts = numshiftbits.(eltype.(barezigs), Ns)
    left_fallback = fallback(zigs[1])
    right_fallback = fallback(zigs[end])
    at = AliasTable(_p)

    CompositeZiggurat(barezigs, pdf, left_fallback, right_fallback, shifts, at)
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
get_subdomains(f, domain) = _get_subdomains(f, promote(float.(domain)...))
get_subdomains(f, domain::NTuple{N,<:AbstractFloat}) where {N} = _get_subdomains(f, domain)
get_subdomains(f, domain::AbstractArray{<:AbstractFloat}) = _get_subdomains(f, domain)
function _get_subdomains(f, domain)
    sd = Vector{Vector{eltype(domain)}}(undef, length(domain)-1)
    for i in eachindex(sd)
        sd[i] = Vector{eltype(domain)}(undef, 2)
    end

    sd[1][1] = domain[1]
    for (i, x) in enumerate(domain[2:(end - 1)])
        if isapprox(f(prevfloat(x)), f(x); atol = sqrt(eps(x)))
            sd[i][2] = x
        else
            sd[i][2] = prevfloat(x)
        end

        if isapprox(f(nextfloat(x)), f(x); atol = sqrt(eps(x)))
            sd[i + 1][1] = x
        else
            sd[i + 1][1] = nextfloat(x)
        end
    end
    sd[end][2] = domain[end]
    sd
end

Random.eltype(::Type{<:CompositeZiggurat{X}}) where {X} = X
Ytype(::Type{<:CompositeZiggurat{X,Y}}) where {X,Y} = Y
Ytype(::CompositeZiggurat{X,Y}) where {X,Y} = Y

function Base.rand(rng::AbstractRNG, zig_sampler::Random.SamplerTrivial{<:CompositeZiggurat})
    @inbounds begin
        zigs = zig_sampler[].zigs
        shifts = zig_sampler[].shifts
        pdf = zig_sampler[].pdf
        left_fallback = zig_sampler[].left_fallback
        right_fallback = zig_sampler[].right_fallback
        at = zig_sampler[].at

        i = rand(rng, at)
        if i == 1
            zigsample(rng, shifts[i], zigs[i], pdf, left_fallback)
        elseif i == length(zigs)
            zigsample(rng, shifts[i], zigs[i], pdf, right_fallback)
        else
            zigsample(rng, shifts[i], zigs[i], pdf, nothing)
        end
    end
end
