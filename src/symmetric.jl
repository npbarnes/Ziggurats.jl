abstract type SymmetricZiggurat{X,Y,LM} <: Ziggurat{X,Y} end

struct BellZiggurat{X,Y,LM,Z<:MonotonicZiggurat{X,Y}} <: SymmetricZiggurat{X,Y,LM}
    mzig::Z

    BellZiggurat(z::MonotonicZiggurat{X,Y}, LM) where {X,Y} = new{X,Y,LM,typeof(z)}(z)

    function BellZiggurat(z::MonotonicZiggurat{X,Y}) where {X,Y}
        LM = layermask_signed(eltype(z), numlayers(z))
        BellZiggurat(z, LM)
    end
end

mask(::SymmetricZiggurat{X,Y,LM}) where {X,Y,LM} = LM

widths(z::SymmetricZiggurat) = widths(z.mzig)
numlayers(z::SymmetricZiggurat) = numlayers(z.mzig)
layerratios(z::SymmetricZiggurat) = layerratios(z.mzig)
heights(z::SymmetricZiggurat) = heights(z.mzig)
highside(z::SymmetricZiggurat) = highside(z.mzig)
lowside(z::SymmetricZiggurat) = lowside(z.mzig)
density(z::SymmetricZiggurat) = density(z.mzig)
fallback(z::SymmetricZiggurat) = fallback(z.mzig)
xarray(z::SymmetricZiggurat) = xarray(z.mzig)

"""
    BellZiggurat(pdf, half_domain, [N]; [ipdf, tailarea, fallback, ...])

Constructs a high performance sampler for a univariate, unimodal, symmetric probability
distributions (a.k.a. bell-shaped distributions) defined by a probability density function,
`pdf`. The pdf must be bell shaped, and must not diverge to inifinity, but may otherwise be
arbitrary - including discontinuous functions. Generate random numbers by passing the
returned ziggurat object to Julia's `rand` or `rand!` functions.

The domain is given as a half-domain, where one endpoint is the mode. Note that the `pdf`
function must be monotonic on the half-domain (a consequence of being unimodal). The `pdf`
function will not be evaluated outside the given half-domain.

The arguments are the same as `monotonic_ziggurat()`. Keep in mind that the
`ipdf`, `tailarea`, and `fallback` arguments are one sided and all must agree on which side
with each other and with the half-domain.
"""
function BellZiggurat(
    pdf,
    half_domain,
    N = nothing;
    ipdf = nothing,
    tailarea = nothing,
    cdf = nothing,
    ccdf = nothing,
    fallback = nothing
)
    z = monotonic_ziggurat(pdf, half_domain, N; ipdf, tailarea, fallback, cdf, ccdf)
    BellZiggurat(z)
end

# The where clause is required to force method specialization
@noinline function bellzigsample_unlikely(rng, w, k, y, mb, am, pdf::F, fb::Nothing, LM, l, x, flip) where {F}
    @inbounds begin
        # Check density, by evaluating the pdf in its proper half-domain
        xx = x_in_domain(x, mb, am, flip)
        yy = (y[l + 1] - y[l]) * rand(rng, eltype(y)) + y[l]

        if yy < pdf(xx)
            return x
        end

        # reject sample and retry
        bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM)
    end
end

# The where clause is required to force method specialization
@noinline function bellzigsample_unlikely(rng, w, k, y, mb, am, pdf::F, fb::FB, LM, l, x, flip) where {F,FB}
    @inbounds begin
        if l == 1
            # unbounded tail fallback
            x2 = w[2] * significand_bitmask(eltype(w)) + mb
            result = convert(eltype(w), fb(rng, x2))
            return ifelse(flip, 2mb-result, result)
        end

        # Check density, by evaluating the pdf in its proper half-domain
        xx = x_in_domain(x, mb, am, flip)
        yy = (y[l + 1] - y[l]) * rand(rng, eltype(y)) + y[l]

        if yy < pdf(xx)
            return x
        end

        # reject sample and retry
        bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM)
    end
end

# The where clause is required to force method specialization
@inline function bellzigsample(rng, w, k, y, mb, am, pdf::F, fb::FB, LM, l, u, s) where {F,FB}
    @inbounds begin
        x = ifelse(s, -u, u) * w[l] + mb
        if u <= k[l]
            return x
        end
        bellzigsample_unlikely(rng, w, k, y, mb, am, pdf, fb, LM, l, x, s)
    end
end

# The where clause is required to force method specialization
@inline function bellzigsample(
    rng,
    w::AbstractArray{<:FloatXX},
    k,
    y,
    mb,
    am,
    pdf::F,
    fb::FB,
    LM::Unsigned
) where {F,FB}
    r = rand(rng, corresponding_uint(eltype(w)))
    bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM, r)
end

# This method is identical to the one above, but it's needed to resolve a method ambiguity
@inline function bellzigsample(rng, w::AbstractArray{<:FloatXX}, k, y, mb, am, pdf::F, fb::FB, LM::Nothing) where {F,FB}
    r = rand(rng, corresponding_uint(eltype(w)))
    bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM, r)
end

# The where clause is required to force method specialization
@inline function bellzigsample(
    rng,
    w::AbstractArray{<:FloatXX},
    k,
    y,
    mb,
    am,
    pdf::F,
    fb::FB,
    LM::Nothing,
    r
) where {F,FB}
    l = rand(rng, 1:(length(w) - 1))
    u = signed(r >>> shiftbits(eltype(w)))
    s = r % Bool
    bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM, l, u, s)
end

# The where clause is required to force method specialization
@inline function bellzigsample(
    rng,
    w::AbstractArray{<:FloatXX},
    k,
    y,
    mb,
    am,
    pdf::F,
    fb::FB,
    LM::Unsigned,
    r
) where {F,FB}
    l = layer_bits_signed(eltype(w), LM, r) + 1
    u = signed(r >>> shiftbits(eltype(w)))
    s = r % Bool
    bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM, l, u, s)
end

# The where clause is required to force method specialization
@inline function bellzigsample(rng, w, k, y, mb, am, pdf::F, fb::FB, LM::Nothing) where {F,FB}
    l = rand(rng, 1:(length(w) - 1))
    u = rand(rng, eltype(w))
    s = rand(rng, Bool)
    bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM, l, u, s)
end

@inline function Base.rand(rng::AbstractRNG, zig_sampler::Random.SamplerTrivial{<:BellZiggurat})
    z = zig_sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    am = lowside(z)
    pdf = density(z)
    fb = fallback(z)
    LM = mask(z)
    bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM)
end

@inline function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:BellZiggurat{X}}
) where {X<:FloatXX}
    z = s[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    am = lowside(z)
    pdf = density(z)
    fb = fallback(z)
    LM = mask(z)
    if length(A) < 7 # TODO: Tune this number
        for i in eachindex(A)
            @inbounds A[i] = rand(rng, s)
        end
    else
        T = corresponding_uint(X)
        GC.@preserve A rand!(rng, Random.UnsafeView{T}(pointer(A), length(A)))

        for i in eachindex(A)
            @inbounds r = reinterpret(T, A[i])
            @inbounds A[i] = bellzigsample(rng, w, k, y, mb, am, pdf, fb, LM, r)
        end
    end
    A
end
