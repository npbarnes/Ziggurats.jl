function between(a, b, x)
    l, r = minmax(a, b)
    l <= x <= r
end

struct Regularized{T}
    a::Vector{T}
end

regularize(domain::Regularized) = domain
function regularize(domain)
    reg = unique(sort!(collect(float.(promote(domain...)))))
    if length(reg) < 2
        error("empty domain. The domain needs at least two distinct points to mark the boundaries.")
    end
    Regularized(reg)
end

Base.getindex(r::Regularized, i) = r.a[i]
Base.extrema(r::Regularized) = extrema(r.a)
Base.length(r::Regularized) = length(r.a)
Base.iterate(r::Regularized) = iterate(r.a)
Base.iterate(r::Regularized, s) = iterate(r.a, s)
Base.firstindex(r::Regularized) = firstindex(r.a)
Base.lastindex(r::Regularized) = lastindex(r.a)
Base.first(r::Regularized) = first(r.a)
Base.last(r::Regularized) = last(r.a)
Base.eltype(r::Regularized) = eltype(r.a)
