struct InverseFunc{X,Y,F}
    f::F
    bigf::Y
    smallf::Y
    smallx::X
    a::X
    b::X

    function InverseFunc{X,Y}(f, (a, b)) where {X,Y}
        InverseFunc{X,Y}(f, a, b)
    end

    """
        InverseFunc{X,Y}(f, a, b)
        InverseFunc{X,Y}(f, (a,b))

    Given a monotonic, non-negative function, `f`, and its domain, `(a,b)`, return a callable
    object that inverts `f` using root finding.

    The result is undefined for any function that is negative or non-monotonic anywhere on the
    given domain. When `0 <= y < min(f(a),f(b))`, `InverseFunc{X,Y}(f, a, b)(y)` returns
    `argmin(f, (a,b))`. This is because increasing functions are treated as being zero for
    `x < a`, and decreasing functions are treated as being zero for `x > b` (a constant
    function is taken to be decreasing). Therefore, the domain of the inverse is
    `0 <= y <= max(f(a),f(b))`. When y is outside that range, an error is raised.

    # Example
    ```julia-repl
    julia> import Ziggurats.InverseFunc

    julia> inverse_func = InverseFunc{Float64,Float64}(x->cos(x) + 2x, (0,3));

    julia> inverse_func(3)
    1.4296716508827783
    ```
    """
    function InverseFunc{X,Y}(f, a, b) where {X,Y}
        a, b = minmax(a, b)
        fa = convert(Y, f(a))
        fb = convert(Y, f(b))
        if isnan(fa) || isnan(fb)
            error("f returns NaN on at least one end point of the domain. The function f must \
            return a number for all values in its domain to make a sensible inverse function.")
        end

        smallx, smallf, bigf = fa < fb ? (a, fa, fb) : (b, fb, fa)

        new{X,Y,typeof(f)}(f, bigf, smallf, smallx, a, b)
    end
end

function (self::InverseFunc{X,Y})(y) where {X,Y}
    if y > self.bigf || y < zero(y)
        throw(ArgumentError("No inverse exists when y > max(f(a), f(b)) or y <= 0."))
    elseif y <= self.smallf
        return self.smallx
    end
    g = let f=self.f, y=y
        x -> convert(Y, f(x)) - y
    end

    find_zero(g, (self.a, self.b))::X
end
