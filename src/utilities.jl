function between(a, b, x)
    l, r = minmax(a, b)
    l <= x <= r
end

struct RollingPairs{T}
    itr::T
    function RollingPairs(itr::T) where {T}
        # Check that the length is zero or 2+.
        result1 = iterate(itr)
        if isnothing(result1)
            return new{T}(itr)
        end

        (_, s) = result1
        result2 = iterate(itr, s)
        if isnothing(result2)
            error("a RollingPairs iterator cannot be formed from an iterator with just one element. There must be zero, 2, or more elements.")
        end

        new{T}(itr)
    end
end

function Base.iterate(r::RollingPairs)
    result1 = iterate(r.itr)
    if result1 === nothing
        return nothing
    end

    val1, state1 = result1

    result2 = iterate(r.itr, state1)
    if result2 === nothing
        error("a RollingPairs iterator cannot be formed from an iterator with just one element. There must be zero, 2, or more elements.")
    end

    val2, state2 = result2

    return ((val1, val2), (val2, state2))
end

function Base.iterate(r::RollingPairs, (prev_val, state))
    result = iterate(r.itr, state)
    if result === nothing
        return nothing
    end

    next_val, next_state = result
    return ((prev_val, next_val), (next_val, next_state))
end

function Base.length(r::RollingPairs)
    len = length(r.itr)
    if len == 0
        return len
    end

    if len == 1
        error("a RollingPairs iterator cannot be formed from an iterator with just one element. There must be zero, 2, or more elements.")
    end

    return len - 1
end

Base.eltype(r::RollingPairs) = Tuple{eltype(r.itr),eltype(r.itr)}

const ManifestlyListOfScalar = Union{
    AbstractVector{<:Number},
    Tuple{Vararg{Number}} # allows empty tuples
}
const ManifestlyListOfListOfScalar = Union{
    AbstractVector{<:ManifestlyListOfScalar},
    Tuple{ManifestlyListOfScalar,Vararg{ManifestlyListOfScalar}} # does not allow empty tuples
}

islistofscalar(::ManifestlyListOfScalar) = true
islistofscalar(::ManifestlyListOfListOfScalar) = false
islistofscalar(a) = _islistofscalar(a)

_islistofscalar(::Number) = false
function _islistofscalar(a)
    if !applicable(iterate, a)
        return false
    end
    for x in a
        if !(x isa Number)
            return false
        end
    end
    return true
end

struct Scalar end
struct ListOfScalar end
struct ListOfListOfScalar end
struct SomethingElse end

classify(::ManifestlyListOfScalar) = ListOfScalar()
classify(::ManifestlyListOfListOfScalar) = ListOfListOfScalar()
classify(a) = _classify(a)

_classify(::Number) = Scalar()
function _classify(a)
    if islistofscalar(a)
        return ListOfScalar()
    elseif applicable(iterate, a) && all(islistofscalar, a)
        return ListOfListOfScalar()
    else
        return SomethingElse()
    end
end

best_float_type(a) = _best_float_type(a)
@generated function best_float_type(a::Tuple)
    types = a.parameters
    if all(x -> isconcretetype(eltype(x)), types)
        T = float(mapreduce(eltype, promote_type, types))
        return :($T)
    end

    :(_best_float_type(a))
end

function _best_float_type(a)
    # works on both ListOfScalar and ListOfListOfScalar because eltype of a Number type is itself
    # and flatten doesn't affect lists that are already flat.
    # i.e. eltype(eltype(::Vector{Float64})) == Float64, and flatten(a::Vector{Float64}) == a
    T = eltype(eltype(a))
    if isconcretetype(T)
        return float(T)
    end

    float(mapreduce(typeof, promote_type, Iterators.flatten(a)))
end

function validate(_, domain)
    error("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got $domain.")
    nothing
end

function validate(::ListOfScalar, domain)
    if !(classify(domain) isa ListOfScalar)
        error("expected a contiguous domain to be an iterator of Number values. Got $domain.")
    end

    if !issorted(domain)
        error("the critical points that make up a contiguous domain must be in increasing order. Got $domain.")
    end

    if !allunique(domain)
        error("duplicate values found in a contiguous domain = $domain.")
    end

    if length(domain) < 2
        error("a contiguous domain needs at least two distinct points to mark the boundaries. Got $domain.")
    end

    nothing
end

function is_overlapping(domain)
    if length(domain) == 1
        return false
    end

    if issorted(domain; by = first)
        d = domain
    else
        d = sort!(collect(domain); by = first)
    end

    for (a, b) in RollingPairs(d)
        if last(a) > first(b)
            return true
        end
    end

    return false
end

function validate(::ListOfListOfScalar, domain)
    foreach(d -> validate(ListOfScalar(), d), domain)

    if is_overlapping(domain)
        error("the contiguous domains that make up `domain` must be non-overlapping. Got $domain.")
    end

    nothing
end

struct CanonicalArray{T<:AbstractFloat} <: AbstractVector{Tuple{T,T}}
    arr::Vector{Tuple{T,T}}
    function CanonicalArray(arr::Vector{Tuple{T,T}}) where {T<:AbstractFloat}
        validate(ListOfListOfScalar(), arr)
        if issorted(arr; by = first)
            s = arr
        else
            s = sort(arr; by = first)
        end
        new{T}(s)
    end
end

Base.getindex(c::CanonicalArray, i) = c.arr[i]
Base.size(c::CanonicalArray) = size(c.arr)
Base.IndexStyle(::Type{<:CanonicalArray}) = IndexLinear()

canonicalize(domain::CanonicalArray) = domain
canonicalize(::Type{A}, domain::CanonicalArray{A}) where {A} = domain

function canonicalize(domain::Vector{Tuple{T,T}}) where {T<:AbstractFloat}
    t = classify(domain)
    if t isa Scalar || t isa SomethingElse
        error("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got $domain.")
    end
    validate(t, domain)
    CanonicalArray(domain)
end
canonicalize(::Type{A}, domain::Vector{Tuple{A,A}}) where {A<:AbstractFloat} = canonicalize(domain)::CanonicalArray{A}

function canonicalize(domain)
    t = classify(domain)
    if t isa Scalar || t isa SomethingElse
        error("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got $domain.")
    end
    validate(t, domain)
    CanonicalArray(_canonicalize(t, domain))
end

function canonicalize(::Type{A}, domain)::CanonicalArray{A} where {A}
    t = classify(domain)
    if t isa Scalar || t isa SomethingElse
        error("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got $domain.")
    end
    validate(t, domain)
    CanonicalArray(_canonicalize(A, t, domain))
end

function _canonicalize(t, domain)
    T = best_float_type(domain)
    _canonicalize(T, t, domain)
end

function _canonicalize(T, ::ListOfScalar, domain)
    result = Vector{Tuple{T,T}}(undef, length(domain) - 1)
    for (i, p) in zip(eachindex(result), RollingPairs(domain))
        result[i] = p
    end

    sort!(result; by = first)
end

function _canonicalize(T, ::ListOfListOfScalar, domain)
    num_subdomains = sum(l -> length(l)-1, domain)
    result = Vector{Tuple{T,T}}(undef, num_subdomains)
    i = 1
    for group in domain
        for p in RollingPairs(group)
            result[i] = p
            i += 1
        end
    end
    @assert i == num_subdomains + 1

    sort!(result; by = first)
end
