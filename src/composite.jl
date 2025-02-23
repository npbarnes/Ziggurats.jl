struct SymmetricZiggurat{X, Z<:MonotonicZiggurat{X}}
    monotonic::Z
    center::X
end

function SymmetricZiggurat(f, domain, N; ipdf_left=nothing, ipdf_right=nothing)
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

    monotonic = monotonic_ziggurat(f, halfdomain, N, ipdf)

    SymmetricZiggurat(monotonic, center)
end

function Base.rand(
    rng::AbstractRNG,
    zig_sampler::Random.SamplerTrivial{SymmetricZiggurat}
)
    monotonic = zig_sampler[].monotonic
    center = zig_sampler[].center

    r = rand(rng, monotonic)
    ifelse(rand(rng, Bool), center + r, center - r)
end

