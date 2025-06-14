"""
    inverse(f, d, args...; kwargs...)

Given a monotonic, non-negative function, `f`, and its domain, `d`, and return a function
that inverts `f` using root finding. 

The result is undefined for any function that is negative or non-monotonic anywhere on the
given domain. When `0 <= y < min(f(a),f(b))`, `inverse(f, (a,b))(y)` returns
`argmin(f, (a,b))`. This is because increasing functions are treated as being zero for
`x < a`, and decreasing functions are treated as being zero for `x > b` (a constant
function is taken to be decreasing). Therefore, the domain of the inverse is
`0 <= y <= max(f(a),f(b))`. When y is outside that range, an error is raised. The domain
may be provided as a pair, (a,b), but it may also include intermediate points, like
`[1,2,3]`. Intermediate points are ignored. Additional arguments are passed along to
`Roots.find_zero`.

# Example
```julia-repl
julia> import ZigguratTools.inverse

julia> inverse(x->cos(x) + 2x, (0,3))(3)
1.4296716508827783
```
"""
function inverse(f, domain, args...; kwargs...)
    ex = extrema(domain)
    a, b = promote(float(ex[1]), float(ex[2]))
    fa, fb = f(a), f(b)
    if isnan(fa) || isnan(fb)
        error("f returns NaN on at least one end point of the domain. The function f must \
        return a number for all values in its domain to make a sensible inverse function.")
    end
    smallx, smallf, bigf = fa < fb ? (a, fa, fb) : (b, fb, fa)

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
