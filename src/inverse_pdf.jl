function inverse(f, domain, args...; kwargs...)
    ex = extrema(domain)
    a, b = promote(float(ex[1]), float(ex[2]))
    fa, fb = f(a), f(b)
    smallx, smallf, bigf = fa < fb ? (a, fa, fb) : (b, fb, fa)

    let bigf=bigf, smallf=smallf, smallx=smallx, f=f, args=args, kwargs=kwargs
        y -> begin
            if y > bigf
                throw(ArgumentError("No inverse exists when y > max(f(a), f(b))."))
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
