"""
    inversepdf(f, d, args...; kwargs...)

Given a monotonic, non-negative function, `f`, and its domain, `d`, and return a function
that inverts `f` using root finding. 

The result is undefined for any function that is negative or non-monotonic anywhere on the
given domain. When `0 <= y < min(f(a),f(b))`, `inversepdf(f, (a,b))(y)` returns
`argmin(f, (a,b))`. This is because increasing functions are treated as being zero for
`x < a`, and decreasing functions are treated as being zero for `x > b` (a constant
function is taken to be decreasing). Therefore, the domain of the inverse is
`0 <= y <= max(f(a),f(b))`. When y is outside that range, an error is raised. The domain
may be provided as a pair, (a,b), but it may also include intermediate points, like
`[1,2,3]`. Intermediate points are ignored. Additional arguments are passed along to
`Roots.find_zero`.

# Example
```julia-repl
julia> import Ziggurats.inversepdf

julia> inversepdf(x->cos(x) + 2x, (0,3))(3)
1.4296716508827783
```
"""
function inversepdf(f, domain, args...; kwargs...)
    ex = extrema(domain)
    a, b = promote(float(ex[1]), float(ex[2]))
    fa, fb = f(a), f(b)
    if isnan(fa) || isnan(fb)
        error("f returns NaN on at least one end point of the domain. The function f must \
        return a number for all values in its domain to make a sensible inverse function.")
    end
    smallx, smallf, bigf = fa < fb ? (a, fa, fb) : (b, fb, fa)

    # TODO: Remove closures
    # TODO: check convergence flag.
    # TODO: Specialize functions: generic_inverse, inverse_monotonic, inverse_nonneg_monotonic,
    # inversepdf, inversecdf
    let bigf=bigf, smallf=smallf, smallx=smallx, f=f, args=args, kwargs=kwargs
        y -> begin
            if y > bigf || y < zero(y)
                throw(ArgumentError("No inverse exists when y > max(f(a), f(b)) or y <= 0."))
            elseif y <= smallf
                return smallx
            end
            g = let f=f, y=y
                x -> f(x) - y
            end
            find_zero(g, (a, b), args...; kwargs...)
        end
    end
end

abstract type AbstractInverse end

"""
    Inverse(f, d)

Invert monotonic function `f` on domain `d`.
"""
struct Inverse{F,X,Y,A,KW} <: AbstractInverse
    f::F
    low::X
    high::X
    f_low::Y
    f_high::Y
    default_args::A
    default_kwargs::KW
end

callfunc(i::AbstractInverse, x) = i.f(x)
argmin(i::AbstractInverse) = i.low
min(i::AbstractInverse) = i.f_low
argmax(i::AbstractInverse) = i.high
max(i::AbstractInverse) = i.f_high
default_args(i::AbstractInverse) = i.default_args
default_kwargs(i::AbstractInverse) = i.default_kwargs
xdomain(i::AbstractInverse) = minmax(argmin(i), argmax(i))
ydomain(i::AbstractInverse) = argmin(i) < argmax(i) ? (min(i), max(i)) : (max(i), min(i))

function inverse_monotonic(f, domain, args...; kwargs...)
    a, b = extrema(regularize(domain))
    fa, fb = f(a), f(b)

    if isnan(fa) || isnan(fb)
        error("f returns NaN on at least one end point of the domain. The function f must \
        return a number for all values in its domain to make a sensible inverse function.")
    end

    # Constant functions are considered decreasing
    low, high, f_low, f_high = fa < fb ? (a, b, fa, fb) : (b, a, fb, fa)

    Inverse(f, low, high, f_low, f_high, args, kwargs)
end

function (inv::Inverse)(y, args...; kwargs...)
    if !between(y, inv.f_low, inv.f_high)
        error("No inverse exists. y must be between f(a) and f(b) for domain (a,b).")
    end
    if isempty(args)
        args = inv.default_args
    end
    if isempty(kwargs)
        kwargs = inv.default_kwargs
    end
    find_zero(x -> inv.f(x) - y, (a, b), args...; kwargs...)
end
