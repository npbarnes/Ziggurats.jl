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
        throw(DomainError(q, "the argument to quantile must be between zero and one inclusive."))
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

function StatsBase.mode(::Union{SteppedExponential{T},Type{SteppedExponential{T}}}) where {T}
    zero(T)
end

Distributions.params(d::SteppedExponential) = (d.θ,)

function ipdf_SteppedExponential(d::SteppedExponential)
    k = rate(d)
    let k=k
        y -> floor(-1/k*log(y/(1-exp(-k)))) + 1
    end
end

mutable struct SteppedExpFallback{X,AT}
    const d::SteppedExponential{X}
    currx::X
    at::AT
    function SteppedExpFallback(d)
        x = zero(eltype(d))
        at = _newaliastable(d, x)
        new{eltype(d),typeof(at)}(d, x, at)
    end
end

function (fb::SteppedExpFallback)(rng, x)
    if x != fb.currx
        fb.currx = x
        fb.at = _newaliastable(fb.d, x)
    end

    bin = rand(rng, fb.at)
    if bin == 1
        if fb.currx == ceil(fb.currx)
            Δ = one(fb.currx)
        else
            Δ = (ceil(fb.currx) - fb.currx)
        end
        return Δ*rand(rng, eltype(fb.d)) + fb.currx
    else
        return floor(fb.currx) + bin - rand(rng, eltype(fb.d))
    end
end

function _newaliastable(d, x)
    start_bin = floor(x)
    start_p = pdf(d, x) * (1 - (x - floor(x)))
    final_bin = ceil(Int, quantile(d, 0.999_999_999))
    probabilities = [start_p; pdf.(d, (start_bin + 1):final_bin)]
    AliasTable(probabilities)
end
