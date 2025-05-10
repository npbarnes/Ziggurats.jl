using Roots

function inverse(f, domain, args...; kwargs...)
    ex = extrema(domain)
    a,b = promote(float(ex[1]), float(ex[2]))
    fa, fb = f(a), f(b)
    smallx, smallf, bigf = fa < fb ? (a, fa, fb) : (b, fb, fa)
    y -> begin
        if y > bigf
            throw(ArgumentError("No inverse exists when y > max(f(a), f(b))."))
        elseif y <= smallf
            return smallx
        end
        find_zero(x -> f(x) - y, (a,b), args...; kwargs...)
    end
end
