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

# find middle of (a,b) with convention that
# * if a, b finite, they are made non-finite
# if a,b of different signs, middle is 0
# middle falls back to a/2 + b/2, but
# for Float64 values, middle is over the
# reinterpreted unsigned integer.
function _middle(x, y)
    a = isinf(x) ? nextfloat(x) : x
    b = isinf(y) ? prevfloat(y) : y
    if sign(a) * sign(b) < 0
        return zero(a)
    else
        __middle(a, b)
    end
end
## find middle assuming a,b same sign, finite
## Alternative "mean" definition that operates on the binary representation
## of a float. Using this definition, bisection will never take more than
## 64 steps (over Float64)
__middle(x::Float64, y::Float64) = __middle(Float64, UInt64, x, y)
__middle(x::Float32, y::Float32) = __middle(Float32, UInt32, x, y)
__middle(x::Float16, y::Float16) = __middle(Float16, UInt16, x, y)

function __middle(T, S, x, y)
    # Use the usual float rules for combining non-finite numbers
    # do division over unsigned integers with bit shift
    xint = reinterpret(S, abs(x))
    yint = reinterpret(S, abs(y))
    mid = (xint + yint) >> 1

    # reinterpret in original floating point
    sign(x + y) * reinterpret(T, mid)
end
