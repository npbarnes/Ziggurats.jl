struct SymmetricZiggurat{X,Z<:MonotonicZiggurat{X}}
    monotonic::Z
    center::X
end

function SymmetricZiggurat(f, domain, N; ipdf_left = nothing, ipdf_right = nothing)
    domian = regularize_domain(domain)

    if isinf(domain[1]) != isinf(domain[2])
        error("invalid domain. Symmetric distributions must have a bounded domain or (-Inf, Inf), got $domain.")
    end

    if isinf(domain[1]) && isinf(domain[2])
        center = zero(domain[1])
    else
        center = (domain[1] + domain[2])/2
    end

    if ipdf_left !== nothing
        halfdomain = (domain[1], center)
        ipdf = ipdf_left
    elseif ipdf_right !== nothing
        halfdomain = (center, domain[2])
        ipdf = ipdf_right
    else
        halfdomain = (center, domain[2])
        ipdf = inverse(f, halfdomain)
    end

    monotonic = monotonic_ziggurat(f, halfdomain, N; ipdf)

    SymmetricZiggurat(monotonic, center)
end

function Base.rand(
    rng::AbstractRNG,
    zig_sampler::Random.SamplerTrivial{<:SymmetricZiggurat}
)
    monotonic = zig_sampler[].monotonic
    center = zig_sampler[].center

    r = rand(rng, monotonic)
    ifelse(rand(rng, Bool), center + r, center - r)
end

struct CompositeZiggurat{Z<:Tuple,AT<:AliasTable}
    zigs::Z
    at::AT
end

"""
    get_subdomains(f, domain)

Convert a list of subdomain boundaries -- `domain` -- into a list of pairs representing subdomains.
Attempts to detect discontinuities in f, and account for them.

# Examples

```julia-repl
julia> ZigguratTools.get_subdomains(cos, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, 0.0]
 [0.0, 1.0]

julia> ZigguratTools.get_subdomains(x -> x>0 ? 1.0 : 0.0, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, 0.0]
 [5.0e-324, 1.0]

julia> ZigguratTools.get_subdomains(x -> x>=0 ? 1.0 : 0.0, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, -5.0e-324]
 [0.0, 1.0]

julia> ZigguratTools.get_subdomains(sign, [-1, 0, 1])
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

function monotonic_segments(f, domain)
    a, b = extrema(regularize_domain(domain))
    df = x -> ForwardDiff.derivative(f, x)

    unique([a; find_zeros(df, a, b); b])
end

# TODO: add option to autodetect monotonic subdomains.
function CompositeZiggurat(pdf, domain, N::Integer; kwargs...)
    Ns = fill(N, length(domain)-1)
    CompositeZiggurat(pdf, domain, Ns; kwargs...)
end

function CompositeZiggurat(
    pdf,
    domain,
    Ns;
    ipdfs = nothing,
    cdf = nothing,
    ccdf = nothing,
    p = nothing
)
    domain = regularize_domain(domain)

    a, b = extrema(domain)
    subdomains = get_subdomains(pdf, domain)

    if length(subdomains) != length(Ns)
        throw(ArgumentError("N must be either an Integer or an iterable with length equal to the number of subdomains."))
    end
    if ipdfs === nothing
        ipdfs = [inverse(pdf, d) for d in subdomains]
    end
    if cdf === nothing
        val = pdf(Roots.__middle(a, b))
        domain_type = typeof(domain[1])
        range_type = typeof(val)
        error_type = typeof(norm(val))
        segbuf = alloc_segbuf(domain_type, range_type, error_type)

        cdf = let pdf=pdf, a=a, segbuf=segbuf
            # Handle the lower boundary explicitly since quadgk will produce NaN.
            x -> if x == a
                zero(domain_type)
            else
                abs(quadgk(pdf, a, max(a, x); segbuf)[1])
            end
        end
    end
    if ccdf === nothing
        val = pdf(Roots.__middle(a, b))
        domain_type = typeof(domain[1])
        range_type = typeof(val)
        error_type = typeof(norm(val))
        segbuf = alloc_segbuf(domain_type, range_type, error_type)

        ccdf = let pdf=pdf, b=b, segbuf=segbuf
            # Handle Inf explicitly since quadgk will produce NaN.
            x -> x == Inf ? zero(domain_type) : abs(quadgk(pdf, min(b, x), b; segbuf)[1])
        end
    end
    if p === nothing
        _p = [cdf(b) - cdf(a) for (a, b) in subdomains]
    else
        _p = similar(p, Float64)
        for i in eachindex(p)
            if p[i] === nothing
                a, b = subdomains[i]
                _p[i] = cdf(b) - cdf(a)
            end
        end
    end

    if !(length(domain) - 1 == length(Ns) == length(ipdfs) == length(_p))
        throw(ArgumentError("Ns, ipdfs, and p must all have a length one less than the length of domain."))
    end

    # TODO: Handle fallback generators
    zig_gen =
        i -> begin
            monotonic_ziggurat(pdf, subdomains[i], Ns[i]; ipdf = ipdfs[i], cdf, ccdf)
        end

    zigs = ntuple(zig_gen, length(subdomains))

    at = AliasTable(_p)

    CompositeZiggurat(zigs, at)
end

function Base.rand(
    rng::AbstractRNG,
    zig_sampler::Random.SamplerTrivial{<:CompositeZiggurat}
)
    zigs = zig_sampler[].zigs
    at = zig_sampler[].at

    i = rand(rng, at)
    rand(rng, zigs[i])
end
