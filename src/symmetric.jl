struct UnimodalZiggurat{X,Y,LM,Z<:MonotonicZiggurat{X,Y}} <: Ziggurat{X,Y}
    mzig::Z

    UnimodalZiggurat(z::MonotonicZiggurat{X,Y}, LM) where {X,Y} = new{X,Y,LM,typeof(z)}(z)

    function UnimodalZiggurat(z::MonotonicZiggurat{X,Y,LM_unsigned}) where {X,Y,LM_unsigned}
        first_bit = typemin(signed(corresponding_uint(X)))
        if LM_unsigned & first_bit == 0
            LM = LM_unsigned << 1
        else
            LM = nothing
        end
        UnimodalZiggurat(z, LM)
    end
end

widths(z::UnimodalZiggurat) = widths(z.mzig)
numlayers(z::UnimodalZiggurat) = numlayers(z.mzig)
layerratios(z::UnimodalZiggurat) = layerratios(z.mzig)
heights(z::UnimodalZiggurat) = heights(z.mzig)
highside(z::UnimodalZiggurat) = highside(z.mzig)
density(z::UnimodalZiggurat) = density(z.mzig)
fallback(z::UnimodalZiggurat) = fallback(z.mzig)
xarray(z::UnimodalZiggurat) = xarray(z.mzig)

@inline function _unimodal_zigsample_floats_masked(rng, w, k, y, mb, pdf::F, fb::FB, LM, r) where {F,FB}
    @inbounds begin
        l = layer_bits_signed(eltype(w), LM, r) + 1
        u = signed(r >>> shiftbits(eltype(w)))
        flip = r % Bool
        x = ifelse(flip, -u, u)*w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(unimodal_zigsample_floats_masked, rng, w, k, y, mb, pdf, fb, LM, l, x, flip)
    end
end

@inline function unimodal_zigsample_floats_masked(rng, w, k, y, mb, pdf::F, fb::FB, LM) where {F,FB}
    r = rand(rng, corresponding_uint(eltype(w)))
    _unimodal_zigsample_floats_masked(rng, w, k, y, mb, pdf, fb, LM, r)
end

@inline function _unimodal_zigsample_floats(rng, w, k, y, mb, pdf::F, fb::FB, LM, r) where {F,FB}
    @inbounds begin
        l = rand(rng, 1:(length(w) - 1))
        u = signed(r >>> shiftbits(eltype(w)))
        flip = r % Bool
        x = ifelse(flip, -u, u)*w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(unimodal_zigsample_floats, rng, w, k, y, mb, pdf, fb, LM, l, x, flip)
    end
end

@inline function unimodal_zigsample_floats(rng, w, k, y, mb, pdf::F, fb::FB, LM) where {F,FB}
    r = rand(rng, corresponding_uint(eltype(w)))
    _unimodal_zigsample_floats(rng, w, k, y, mb, pdf, fb, LM, r)
end

@inline function unimodal_zigsample_general(rng, w, k, y, mb, pdf::F, fb::FB, LM, r) where {F,FB}
    @inbounds begin
        l = rand(rng, 1:(length(w) - 1))
        u = rand(rng, eltype(w))
        flip = r % Bool
        x = ifelse(flip, -u, u)*w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(unimodal_zigsample_general, rng, w, k, y, mb, pdf, fb, LM, l, x, flip)
    end
end

@inline function Base.rand(
    rng::AbstractRNG,
    sampler::Random.SamplerTrivial{<:UnimodalZiggurat{X,Y,LM}}
) where {X<:FloatXX,Y,LM}
    z = sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    unimodal_zigsample_floats_masked(rng, w, k, y, mb, pdf, fb, LM)
end

@inline function Base.rand(
    rng::AbstractRNG,
    sampler::Random.SamplerTrivial{<:UnimodalZiggurat{X,Y,nothing}}
) where {X<:FloatXX,Y}
    z = sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    unimodal_zigsample_floats(rng, w, k, y, mb, pdf, fb, nothing)
end

@inline function Base.rand(rng::AbstractRNG, sampler::Random.SamplerTrivial{<:UnimodalZiggurat{X,Y}}) where {X,Y}
    z = sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    unimodal_zigsample_general(rng, w, k, y, mb, pdf, fb, nothing)
end
