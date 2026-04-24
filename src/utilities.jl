function between(a, b, x)
    l, r = minmax(a, b)
    l <= x <= r
end

struct Regularized{T}
    a::T
end

promote_to_float(a) = float.(promote(a...))
@generated function promote_to_float(a::AbstractVector{X}) where {X}
    if isconcretetype(X) && X <: AbstractFloat
        return :(copy(a))
    elseif isconcretetype(X) && !(X <: AbstractFloat)
        return :(float.(a))
    elseif !isconcretetype(X) && X <: AbstractFloat
        return :(collect(promote(a...)))
    elseif !isconcretetype(X) && !(X <: AbstractFloat)
        if X <: Real
            return :(Vector{Float64}(a))
        else
            return :(collect(float.(promote(a...))))
        end
    else
        error("Unreachable error")
    end
end

# Workaround for no sort(::NTuple) method in the LTS version Julia 1.10.4
sort_maybeinplace(x::Tuple) = Tuple(sort(SVector(x)))
sort_maybeinplace(x::AbstractArray) = sort!(x)

regularize(domain::Regularized) = domain
function regularize(domain)
    reg = domain |> promote_to_float |> sort_maybeinplace
    if !allunique(reg)
        error("duplicate values found in domain.")
    end
    if length(reg) < 2
        error("The domain needs at least two distinct points to mark the boundaries.")
    end
    Regularized(reg)
end

regularize(::Type{A}, domain::Regularized{B}) where {A,B<:AbstractArray{A}} = domain
regularize(::Type{A}, domain::Regularized{B}) where {A,N,B<:NTuple{N,A}} = domain
regularize(a::Type{A}, domain::Regularized) where {A} = throw(MethodError(regularize, (a, domain)))
function regularize(T, domain)
    reg = convert.(T, domain) |> sort_maybeinplace
    if !allunique(reg)
        error("duplicate values found in domain.")
    end
    if length(reg) < 2
        error("The domain needs at least two distinct points to mark the boundaries.")
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
