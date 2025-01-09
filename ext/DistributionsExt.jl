module DistributionsExt

using ZigguratTools
using Distributions

function _ipdf(dist, y)
    ZigguratTools.inverse(Base.Fix1(pdf, dist), extrema(dist), y)
end

function _ipdf(dist::Exponential, y)
    θ = scale(dist)
    (log(θ) - log(y)) / θ
end

# TODO: Add more _ipdf methods for different distributions

function _modaldensity(dist)
    pdf(dist, mode(dist))
end

function _modaldensity(dist::Exponential)
    rate(dist)
end

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

struct TailFallback{T<:Truncated}
    truncdist::T
end

function left_fallback(dist, x)
    TailFallback(truncated(dist, upper=x))
end

function right_fallback(dist, x)
    TailFallback(truncated(dist, lower=x))
end

function (tf::TailFallback)(rng)
    rand(rng, tf.truncdist)
end

ZigguratTools.monotonic_ziggurat(dist::Distribution, N::Integer=256; kwargs...) = monotonic_ziggurat(dist, extrema(dist), N; kwargs...)
function ZigguratTools.monotonic_ziggurat(dist::Distribution, domain, N=256;
    pdf = Base.Fix1(Distributions.pdf, dist),
    ipdf = Base.Fix1(DistributionsExt.ipdf, dist),
    tailarea = nothing,
    fallback_generator = nothing
    )

    if tailarea === nothing
        if isinf(domain[1])
            tailarea = Base.Fix1(cdf, dist)
        elseif isinf(domain[2])
            tailarea = Base.Fix1(ccdf, dist)
        end
    end

    if fallback_generator === nothing
        if isinf(domain[1])
            fallback_generator = Base.Fix1(left_fallback, dist)
        elseif isinf(domain[2])
            fallback_generator = Base.Fix1(right_fallback, dist)
        end
    end

    monotonic_ziggurat(pdf, domain, N; ipdf, tailarea, fallback_generator)
end

end # module