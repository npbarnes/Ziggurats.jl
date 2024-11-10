"""
    inverse(f, domain)

Returns a function that will compute the inverse of f on the domain. Similar to
`y -> inverse(f, domain, y)`.
"""
function inverse(f, (a,b))
    let f = f, domain = promote(float(a), float(b))
        y -> inverse(f, domain, y)
    end
end

"""
    inverse(f, domain, y)

Find the generalized inverse of the non-constant monotonic function f at y in
the domain. If f is decreasing, then the largest x such that f(x) >= y is
returned. If f is increasing then the smallest x such that f(x) >= y is
returned. If there is no solution (i.e. when y is not between f(a) and f(b))
then an error is thrown. If f is constant, then an error is thrown. If
non-monotonicity is detected in f then an error is thrown, but in general the
result is undefined if f is not monotonic. Note that while f cannot be constant
on the whole domain, it need not be strictly monotonic.
"""
inverse(f, (a, b), y) = inverse(f, promote(float(a), float(b)), y)

function inverse(f, domain::NTuple{2,<:AbstractFloat}, y)
    a, b = domain

    if a > b
        error("domain must be an ordered tuple.")
    end

    fa = f(a)
    fb = f(b)

    if fa == fb
        error("f must be non-constant.")
    end

    if y > max(fa, fb)
        error("no solutions.")
    end

    if fa >= fb
        if fb >= y
            b
        else
            _decreasing_inverse(f, fa, fb, a, b, y)
        end
    else
        if fa >= y
            a
        else
            _increasing_inverse(f, fa, fb, a, b, y)
        end
    end
end

"Largest x such that f(x) >= y"
function _decreasing_inverse(f, fa, fb, a, b, y)
    while nextfloat(a) != b
        c = _middle(a, b)
        fc = f(c)

        if !(fb <= fc <= fa)
            error("f must be monotonic.")
        end

        if fc >= y
            a = c
            fa = fc
        else
            b = c
            fb = fc
        end
    end

    return a
end

"Smallest x such that f(x) >= y"
function _increasing_inverse(f, fa, fb, a, b, y)
    while nextfloat(a) != b
        c = _middle(a, b)
        fc = f(c)

        if !(fa <= fc <= fb)
            error("f must be monotonic.")
        end

        if fc >= y
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
    end

    return b
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
