_meanoffset(y, σ2) = √(-2σ2 * log(√(2π * σ2) * y))
ipdf_right(dist::Normal, y) = mean(dist) + _meanoffset(y, var(dist))
ipdf_left(dist::Normal, y) = mean(dist) - _meanoffset(y, var(dist))

function ipdf_right(dist::Truncated{<:Normal}, y)
    ut = dist.untruncated
    mode_pd = pdf(dist, mode(dist))
    if y > mode_pd
        throw(DomainError(
            y,
            "y must be less than or equal to the density at the mode = $mode_pd"
        ))
    end

    y < pdf(dist, maximum(dist)) ? maximum(dist) : ipdf_right(ut, y * dist.tp)
end

function ipdf_left(dist::Truncated{<:Normal}, y)
    ut = dist.untruncated
    mode_pd = pdf(dist, mode(dist))
    if y > mode_pd
        throw(DomainError(
            y,
            "y must be less than or equal to the density at the mode = $mode_pd"
        ))
    end

    y < pdf(dist, minimum(dist)) ? minimum(dist) : ipdf_left(ut, y * dist.tp)
end

ipdf_right(dist::Exponential, y) = -scale(dist) * log(scale(dist) * y)

"""
    max_root(f, domain)

Returns largest x in the domain (inclusive) such that f(x) >= 0. The function f
must be monotonic, and f(a) and f(b) can not both be negative.
"""
function max_root(f, domain::Tuple{Float64,Float64})
    # This is basically a bisection search. The function f might not be
    # continuous or strictly monotonic, but in floating point there will always
    # be a unique largest x where f(x) >= 0. I copied the method of
    # reinterpreting floats as unsigned ints from Roots.jl.

    a, b = domain

    if a > b
        error("domain must be an ordered tuple.")
    end

    fa = f(a)
    fb = f(b)

    if fa < 0 && fb < 0
        error("no solutions, f(a) and f(b) are both negative.")
    end

    if fb >= 0
        return b
    end

    while nextfloat(a) != b
        c = _middle(a,b)
        fc = f(c)

        if fc >= 0
            a = c
            fa = fc
        else
            b = c
            fb = fc
        end
    end

    return a
end

# Reinterpreting floats as unsigned ints avoids some issues with floating point.
# Inspired by Roots.jl
function _middle(x, y)
    if sign(x) * sign(y) < 0
        zero(x)
    else
        __middle(x, y)
    end
end

__middle(x::Float64, y::Float64) = __middle(Float64, UInt64, x, y)
__middle(x::Float32, y::Float32) = __middle(Float32, UInt32, x, y)
__middle(x::Float16, y::Float16) = __middle(Float16, UInt16, x, y)
## fallback for non FloatNN number types
__middle(x::Number, y::Number) = x / 2 + y / 2

function __middle(T, S, x, y)
    xint = reinterpret(S, abs(x))
    yint = reinterpret(S, abs(y))
    mid = (xint + yint) >> 1

    # reinterpret in original floating point
    sign(x + y) * reinterpret(T, mid)
end
