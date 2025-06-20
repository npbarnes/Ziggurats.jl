module DistributionsExt

using ZigguratTools
using Distributions

function _ipdf(dist, y)
    ZigguratTools.inversepdf(Base.Fix1(pdf, dist), extrema(dist))(y)
end

function _ipdf(dist::Exponential, y)
    θ = scale(dist)
    (log(θ) - log(y)) / θ
end

# TODO: Add more _ipdf methods for different distributions

function ipdf(dist::ContinuousUnivariateDistribution, y)
    # assume that the pdf is monotonic and non-constant
    modaldensity = max(pdf(dist, minimum(dist)), pdf(dist, maximum(dist)))

    if y < 0 || y > modaldensity
        throw(DomainError(y, "y is outside the range of the pdf of dist."))
    end

    if y == 0
        # since d is monotonic, it's mode is on the boundary
        if mode(dist) == minimum(dist)
            return maximum(dist)
        else # mode(dist) == maximum(dist)
            return minimum(dist)
        end
    end

    _ipdf(dist, y)
end

struct LeftFallback{T}
    dist::T
end

struct RightFallback{T}
    dist::T
end

function (f::LeftFallback)(rng, x)
    rand(rng, truncated(f.dist; upper = x))
end

function (f::RightFallback)(rng, x)
    rand(rng, truncated(f.dist; lower = x))
end

function ZigguratTools.monotonic_ziggurat(dist::Distribution, N::Integer = 256; kwargs...)
    monotonic_ziggurat(dist, extrema(dist), N; kwargs...)
end

function ZigguratTools.monotonic_ziggurat(
    dist::Distribution,
    domain,
    N = 256;
    pdf = Base.Fix1(Distributions.pdf, dist),
    ipdf = Base.Fix1(DistributionsExt.ipdf, dist),
    tailarea = nothing,
    fallback = nothing
)
    if tailarea === nothing
        if isinf(domain[1])
            tailarea = Base.Fix1(cdf, dist)
        elseif isinf(domain[2])
            tailarea = Base.Fix1(ccdf, dist)
        end
    end

    if fallback === nothing
        if isinf(domain[1])
            fallback = LeftFallback(dist)
        elseif isinf(domain[2])
            fallback = RightFallback(dist)
        end
    end

    monotonic_ziggurat(pdf, domain, N; ipdf, tailarea, fallback)
end

end # module
