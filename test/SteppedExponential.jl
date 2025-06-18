# An unusual distribution to test discontinuities, and non-strict monotonicity.
struct SteppedExponential{T} <: ContinuousUnivariateDistribution
    θ::T
end
SteppedExponential() = SteppedExponential(1.0)

Base.eltype(::SteppedExponential{T}) where {T} = T

Distributions.scale(d::SteppedExponential) = d.θ
Distributions.rate(d::SteppedExponential) = 1 / d.θ

function Distributions.pdf(d::SteppedExponential, x::Real)
    if x < zero(x)
        return zero(eltype(x))
    end

    k = rate(d)
    (1 - exp(-k)) * exp(-k * floor(x))
end

function Distributions.logpdf(d::SteppedExponential, x::Real)
    if x < zero(x)
        return typemin(eltype(d))
    end

    k = rate(d)
    log(1 - exp(-k)) - k * floor(x)
end

function Distributions.cdf(d::SteppedExponential, x::Real)
    if x <= zero(x)
        return zero(eltype(d))
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
        return typemax(eltype(d))
    end
    θ = scale(d)
    k = rate(d)
    m = floor(-θ * log(1 - q))
    (q + exp(-k * m) - 1) / ((1 - exp(-k)) * exp(-k * m)) + m
end

Base.minimum(::Union{SteppedExponential{T},Type{SteppedExponential{T}}}) where {T} = zero(T)
function Base.maximum(::Union{SteppedExponential{T},Type{SteppedExponential{T}}}) where {T}
    typemax(T)
end

Distributions.insupport(::SteppedExponential, x::Real) = x >= 0

function StatsBase.mode(::Union{
    SteppedExponential{T},
    Type{SteppedExponential{T}}
}) where {T}
    zero(T)
end

Distributions.params(d::SteppedExponential) = (d.θ,)
