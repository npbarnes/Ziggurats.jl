abstract type SymmetricZiggurat{X,Y,LM} <: Ziggurat{X,Y} end

struct BellZiggurat{X,Y,LM,Z<:MonotonicZiggurat{X,Y}} <: SymmetricZiggurat{X,Y,LM}
    mzig::Z

    BellZiggurat(z::MonotonicZiggurat{X,Y}, LM) where {X,Y} = new{X,Y,LM,typeof(z)}(z)

    function BellZiggurat(z::MonotonicZiggurat{X,Y}) where {X,Y}
        LM = layermask_signed(eltype(z), numlayers(z))
        BellZiggurat(z, LM)
    end
end

widths(z::SymmetricZiggurat) = widths(z.mzig)
numlayers(z::SymmetricZiggurat) = numlayers(z.mzig)
layerratios(z::SymmetricZiggurat) = layerratios(z.mzig)
heights(z::SymmetricZiggurat) = heights(z.mzig)
highside(z::SymmetricZiggurat) = highside(z.mzig)
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

@inline function _bellzigsample_floats_masked(rng, w, k, y, mb, pdf::F, fb::FB, LM, r) where {F,FB}
    @inbounds begin
        l = layer_bits_signed(eltype(w), LM, r) + 1
        u = signed(r >>> shiftbits(eltype(w)))
        flip = r % Bool
        x = ifelse(flip, -u, u)*w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(bellzigsample_floats_masked, rng, w, k, y, mb, pdf, fb, LM, l, x, flip)
    end
end

@inline function bellzigsample_floats_masked(rng, w, k, y, mb, pdf::F, fb::FB, LM) where {F,FB}
    r = rand(rng, corresponding_uint(eltype(w)))
    _bellzigsample_floats_masked(rng, w, k, y, mb, pdf, fb, LM, r)
end

@inline function _bellzigsample_floats(rng, w, k, y, mb, pdf::F, fb::FB, LM, r) where {F,FB}
    @inbounds begin
        l = rand(rng, 1:(length(w) - 1))
        u = signed(r >>> shiftbits(eltype(w)))
        flip = r % Bool
        x = ifelse(flip, -u, u)*w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(bellzigsample_floats, rng, w, k, y, mb, pdf, fb, LM, l, x, flip)
    end
end

@inline function bellzigsample_floats(rng, w, k, y, mb, pdf::F, fb::FB, LM) where {F,FB}
    r = rand(rng, corresponding_uint(eltype(w)))
    _bellzigsample_floats(rng, w, k, y, mb, pdf, fb, LM, r)
end

@inline function bellzigsample_general(rng, w, k, y, mb, pdf::F, fb::FB, LM) where {F,FB}
    @inbounds begin
        l = rand(rng, 1:(length(w) - 1))
        u = rand(rng, eltype(w))
        flip = rand(rng, Bool)
        x = ifelse(flip, -u, u)*w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(bellzigsample_general, rng, w, k, y, mb, pdf, fb, LM, l, x, flip)
    end
end

@inline function Base.rand(
    rng::AbstractRNG,
    sampler::Random.SamplerTrivial{<:BellZiggurat{X,Y,LM}}
) where {X<:FloatXX,Y,LM}
    z = sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    bellzigsample_floats_masked(rng, w, k, y, mb, pdf, fb, LM)
end

@inline function Base.rand(
    rng::AbstractRNG,
    sampler::Random.SamplerTrivial{<:BellZiggurat{X,Y,nothing}}
) where {X<:FloatXX,Y}
    z = sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    bellzigsample_floats(rng, w, k, y, mb, pdf, fb, nothing)
end

@inline function Base.rand(rng::AbstractRNG, sampler::Random.SamplerTrivial{<:BellZiggurat{X,Y}}) where {X,Y}
    z = sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    bellzigsample_general(rng, w, k, y, mb, pdf, fb, nothing)
end

function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:BellZiggurat{X,Y,LM}}
) where {X<:FloatXX,Y,LM}
    z = s[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    if length(A) < 7 # TODO: Tune this number
        for i in eachindex(A)
            @inbounds A[i] = rand(rng, s)
        end
    else
        T = corresponding_uint(X)
        # UnsafeView is an internal implementation detail of Random.jl
        GC.@preserve A rand!(rng, Random.UnsafeView{T}(pointer(A), length(A)))

        for i in eachindex(A)
            @inbounds r = reinterpret(T, A[i])
            @inbounds A[i] = _bellzigsample_floats_masked(rng, w, k, y, mb, pdf, fb, LM, r)
        end
    end
    A
end
