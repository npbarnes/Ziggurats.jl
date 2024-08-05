# An unusual distribution to test discontinuities, and non-strict monotonicity.
struct SteppedExponential <: ContinuousUnivariateDistribution
    θ::Float64
end
SteppedExponential() = SteppedExponential(1.0)

Distributions.scale(d::SteppedExponential) = d.θ
Distributions.rate(d::SteppedExponential) = 1 / d.θ

function Distributions.pdf(d::SteppedExponential, x::Real)
    if x < zero(x)
        return 0.0
    end

    k = rate(d)
    (1 - exp(-k)) * exp(-k * floor(x))
end

function Distributions.logpdf(d::SteppedExponential, x::Real)
    if x < zero(x)
        return -Inf
    end

    k = rate(d)
    log(1 - exp(-k)) - k * floor(x)
end

function Distributions.cdf(d::SteppedExponential, x::Real)
    if x <= zero(x)
        return 0.0
    end

    fpart, floorx = modf(x)

    k = rate(d)
    ekx = exp(-k * floorx)
    (1 - ekx) + (1 - exp(-k)) * ekx * fpart
end

function Distributions.quantile(d::SteppedExponential, q::Real)
    if !(zero(q) <= q <= oneunit(q))
        throw(DomainError(
            q,
            "the argument to quantile must be between zero and one inclusive."
        ))
    end

    if q == oneunit(q)
        return Inf
    end
    θ = scale(d)
    k = rate(d)
    m = floor(-θ * log(1 - q))
    (q + exp(-k * m) - 1) / ((1 - exp(-k)) * exp(-k * m)) + m
end

Base.minimum(::Union{SteppedExponential,Type{SteppedExponential}}) = 0.0
Base.maximum(::Union{SteppedExponential,Type{SteppedExponential}}) = Inf

Distributions.insupport(::SteppedExponential, x::Real) = x >= 0.0

StatsBase.mode(::Union{SteppedExponential,Type{SteppedExponential}}) = 0.0

Distributions.params(d::SteppedExponential) = (d.θ,)

function ZigguratTools.ipdf_right(d::SteppedExponential, y)
    if y > pdf(d, 0)
        throw(DomainError(y, "y must be between zero and the maximum of the pdf."))
    end

    ceil(scale(d) * (log(pdf(d, 0)) - log(y)))
end
