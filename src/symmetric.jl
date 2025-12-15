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
