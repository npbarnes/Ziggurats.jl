function between(a, b, x)
    l, r = minmax(a, b)
    l <= x <= r
end

struct Regularized{T}
    a::Vector{T}
end

regularize(domain::Regularized) = domain
function regularize(domain)
    Regularized(unique(sort!(collect(promote(float.(domain)...)))))
end

Base.getindex(r::Regularized, i) = r.a[i]
Base.extrema(r::Regularized) = extrema(r.a)
Base.length(r::Regularized) = length(r.a)
Base.iterate(r::Regularized) = iterate(r.a)
Base.iterate(r::Regularized, s) = iterate(r.a, s)
