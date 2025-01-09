"""
    inverse(f, domain)

Returns a function that will compute the inverse of f on the domain. Similar to
`y -> inverse(f, domain, y)`.
"""
function inverse(f, domain; xatol = zero(domain[1]))
    a, b = domain
    let f = f, (a, b, xatol) = promote(float(a), float(b), float(xatol))
        y -> inverse(f, (a, b), y; xatol)
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
function inverse(f, domain, y; xatol = zero(domain[1]))
    a, b = domain
    a, b, xatol = promote(float(a), float(b), float(xatol))
    inverse(f, (a, b), y; xatol)
end

function inverse(
    f,
    domain::NTuple{2,T},
    y;
    xatol::T = zero(domain[1])
) where {T<:AbstractFloat}
    a, b = domain

    if a > b
        error("domain must be an ordered tuple, got domain=$domain.")
    end

    fa = f(a)
    fb = f(b)

    if fa == fb
        error("f must be non-constant, got f(a) == f(b) == $fa.")
    end

    if y > max(fa, fb)
        error("no solutions exist when y > max(f(a), f(b)), got y=$y, f(a)=$fa, and f(b)=$fb.")
    end

    if fa >= fb
        if fb >= y
            b
        else
            _decreasing_inverse(f, fa, fb, a, b, y, xatol)
        end
    else
        if fa >= y
            a
        else
            _increasing_inverse(f, fa, fb, a, b, y, xatol)
        end
    end
end

"Largest x such that f(x) >= y"
function _decreasing_inverse(f, fa, fb, a, b, y, xatol)
    while nextfloat(a) != b && b - a >= xatol
        c = _middle(a, b)
        fc = f(c)

        if !(fb <= fc <= fa)
            error("f must be monotonic, got f($a) = $fa, f($c) = $fc, and f($b) = $fb.")
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
function _increasing_inverse(f, fa, fb, a, b, y, xatol)
    while nextfloat(a) != b && b - a >= xatol
        c = _middle(a, b)
        fc = f(c)

        if !(fa <= fc <= fb)
            error("f must be monotonic, got f($a) = $fa, f($c) = $fc, and f($b) = $fb.")
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
