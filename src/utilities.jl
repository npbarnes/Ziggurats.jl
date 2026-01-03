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
Base.eltype(r::Regularized) = eltype(r.a)

### UnsafeView
# internal array-like type to circumvent the lack of flexibility with reinterpret
# copied from Random.jl, to avoid using stdlib internals
struct UnsafeView{T} <: DenseArray{T,1}
    ptr::Ptr{T}
    len::Int
end

Base.getindex(a::UnsafeView, i::Int) = unsafe_load(a.ptr, i)
Base.setindex!(a::UnsafeView, x, i::Int) = unsafe_store!(a.ptr, x, i)
Base.pointer(a::UnsafeView) = a.ptr
Base.size(a::UnsafeView) = (a.len,)
Base.elsize(::Type{UnsafeView{T}}) where {T} = sizeof(T)
