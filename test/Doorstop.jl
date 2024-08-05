struct Doorstop <: ContinuousUnivariateDistribution
    a::Float64
    b::Float64
    c::Float64
    h::Float64
    p::Float64
    function Doorstop(backside, corner, tip)
        if !(backside <= corner <= tip)
            throw(ArgumentError("backside, corner, and tip must be in increasing order."))
        end
        a, b, c = backside, corner, tip
        h = 1 / ((b - a) + 1 / 2 * (c - b))
        p = h * (b - a)
        new(a, b, c, h, p)
    end
end

function Distributions.pdf(d::Doorstop, x::Real)
    a, b, c, h = d.a, d.b, d.c, d.h
    if a <= x <= b
        h
    elseif b <= x <= c
        h * (1 - (x - b) / (c - b))
    elseif x < a || x > c
        zero(h)
    elseif isnan(x)
        NaN
    else
        error("Unreachable.")
    end
end

Distributions.logpdf(d::Doorstop, x::Real) = log(pdf(d, x))

function Distributions.cdf(d::Doorstop, x::Real)
    a, b, c, h = d.a, d.b, d.c, d.h

    if x <= a
        zero(x)
    elseif a <= x <= b
        h * (x - a)
    elseif b <= x <= c
        h * (x - a - (x - b)^2 / (2 * (c - b)))
    elseif x >= c
        one(x)
    elseif isnan(x)
        NaN
    else
        error("Unreachable.")
    end
end

function Distributions.quantile(d::Doorstop, q::Real)
    a, b, c, h = d.a, d.b, d.c, d.h

    if q < h * (b - a)
        q / h + a
    else
        (c - b) + b - (c - b) * √(1 + 2 * (b - a) / (c - b) - 2 * q / (h * (c - b)))
    end
end

Base.minimum(d::Doorstop) = d.a
Base.maximum(d::Doorstop) = d.c
Distributions.insupport(d::Doorstop, x::Real) = d.a <= x <= d.c

Distributions.params(d::Doorstop) = (d.a, d.b, d.c)
