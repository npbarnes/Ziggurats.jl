struct CompositeZiggurat{
    X,
    Y,
    Ns,
    M,
    Mm2,
    F,
    ZL<:MonotonicZiggurat{X,Y,N,K,NP1,LM,F} where {N,K,NP1,LM},
    ZM<:NTuple{Mm2,BoundedZiggurat{X,Y,N,K,NP1,LM,F} where {N,K,NP1,LM}},
    ZR<:MonotonicZiggurat{X,Y,N,K,NP1,LM,F} where {N,K,NP1,LM},
    AT
} <: Ziggurat{X,Y}
    zl::ZL
    zm::ZM
    zr::ZR
    at::AT

    function CompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT}(
        zl::ZL,
        zm::ZM,
        zr::ZR,
        at
    ) where {X,Y,Ns<:Tuple,M,Mm2,F,ZL,ZM,ZR,AT}
        ns = Ns.parameters
        if numlayers(zl) != ns[1]
            error("the number of layers in one of the ziggurats does not match the number provided in the type parameter.")
        end
        for (z, n) in zip(zm, ns[2:(end - 1)])
            if numlayers(z) != n
                error("the number of layers in one of the ziggurats does not match the number provided in the type parameter.")
            end
        end
        if numlayers(zr) != ns[end]
            error("the number of layers in one of the ziggurats does not match the number provided in the type parameter.")
        end

        if !isinteger(M) || M < 2
            error("the number of subdomains, M, must be an integer greater than one.")
        end

        if Mm2 != M - 2
            error("the type parameter `Mm2` must be exactly `M - 2`.")
        end

        new{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT}(zl, zm, zr, at)
    end
end

density(z::CompositeZiggurat) = density(z.zigs[1])

const BoundedCompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT} =
    CompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT} where {X,Y,Ns,M,Mm2,F,ZL<:BoundedZiggurat,ZM,ZR<:BoundedZiggurat,AT}

const LeftTailCompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT} =
    CompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT} where {X,Y,Ns,M,Mm2,F,ZL<:UnboundedZiggurat,ZM,ZR<:BoundedZiggurat,AT}

const RightTailCompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT} =
    CompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT} where {X,Y,Ns,M,Mm2,F,ZL<:BoundedZiggurat,ZM,ZR<:UnboundedZiggurat,AT}

const TwoTailCompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR,AT} = CompositeZiggurat{
    X,
    Y,
    Ns,
    M,
    Mm2,
    F,
    ZL,
    ZM,
    ZR,
    AT
} where {X,Y,Ns,M,Mm2,F,ZL<:UnboundedZiggurat,ZM,ZR<:UnboundedZiggurat,AT}

tupletypelength(::Type{T}) where {T<:Tuple} = fieldcount(T)

tupletypefirst(::Type{<:Tuple{First,Vararg}}) where {First} = First
tupletypelast(::Type{T}) where {T<:Tuple{Any,Vararg}} = fieldtype(T, fieldcount(T))

for CZ in
    (:BoundedCompositeZiggurat, :LeftTailCompositeZiggurat, :RightTailCompositeZiggurat, :TwoTailCompositeZiggurat)
    @eval begin
        function $(CZ){X,Y}(pdf, domain::Union{NTuple{Mp1,Any},StaticVector{Mp1}}; kwargs...) where {X,Y,Mp1}
            M = Mp1 - 1
            N = default_numlayers(X, Y)
            $(CZ){X,Y,N,M}(pdf, domain; kwargs...)
        end

        function $(CZ){X,Y,Ns}(pdf, domain; kwargs...) where {X,Y,M,Ns<:NTuple{M,Any}}
            $(CZ){X,Y,Ns,M}(pdf, domain; kwargs...)
        end

        function $(CZ){X,Y,nothing,M}(pdf, domain; kwargs...) where {X,Y,M}
            N = default_numlayers(X, Y)
            $(CZ){X,Y,N,M}(pdf, domain; kwargs...)
        end

        function $(CZ){X,Y,Ns,M}(pdf, domain; kwargs...) where {X,Y,Ns<:Tuple,M}
            # Notice that the case `where {X,Y,M,Ns<:Tuple{Vararg{Any,M}}} is defined below.
            # This case is an error when the length of Ns does not equal M.
            error("the Tuple of layer counts, Ns, must have length equal to M.")
        end

        function $(CZ){X,Y,N,M}(pdf, domain; kwargs...) where {X,Y,N,M}
            # This method is called when four type parameters are given, and N is not a
            # subtype of Tuple. It should be an Int value.
            Ns = NTuple{M,N}
            $(CZ){X,Y,Ns,M}(pdf, domain; kwargs...)
        end
    end
end

@generated function tupletype_default_numlayers(X, Y, Ns)
    _X = extracttype(X)
    _Y = extracttype(Y)
    _Ns = extracttype(Ns)

    vals = default_numlayers.(_X, _Y, _Ns.parameters)

    :(Tuple{$(vals)...})
end

function BoundedCompositeZiggurat{X,Y,Ns,M}(
    pdf,
    domain;
    area = nothing,
    cdf = nothing,
    ccdf = nothing,
    ipdfs = nothing,
    p = nothing
) where {X,Y,M,Ns<:Tuple{Vararg{Any,M}}}
    _Ns = tupletype_default_numlayers(X, Y, Ns)
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, _Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)

    zig_L = BoundedZiggurat{X,Y,tupletypefirst(_Ns)}(pdf, subdomains[1, :]; ipdf = _ipdfs[1])
    zig_R = BoundedZiggurat{X,Y,tupletypelast(_Ns)}(pdf, subdomains[end, :]; ipdf = _ipdfs[end])

    Mm2 = M - 2
    F = typeof(density(zig_L))
    ZL = typeof(zig_L)
    ZM = typeof(zigs_middle)
    ZR = typeof(zig_R)
    AT = typeof(at)

    CompositeZiggurat{X,Y,_Ns,M,Mm2,F,ZL,ZM,ZR,AT}(zig_L, zigs_middle, zig_R, at)
end

function LeftTailCompositeZiggurat{X,Y,Ns,M}(
    pdf,
    domain;
    ipdfs = nothing,
    area = nothing,
    cdf = nothing,
    ccdf = nothing,
    left_fallback = nothing,
    p = nothing
) where {X,Y,M,Ns<:Tuple{Vararg{Any,M}}}
    _Ns = tupletype_default_numlayers(X, Y, Ns)
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, _Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)

    zig_L = UnboundedZiggurat{X,Y,tupletypefirst(_Ns)}(
        pdf,
        subdomains[1, :];
        ipdf = _ipdfs[1],
        tailarea = cdf,
        fallback = left_fallback
    )
    zig_R = BoundedZiggurat{X,Y,tupletypelast(_Ns)}(pdf, subdomains[end, :]; ipdf = _ipdfs[end])

    Mm2 = M - 2
    F = typeof(density(zig_L))
    ZL = typeof(zig_L)
    ZM = typeof(zigs_middle)
    ZR = typeof(zig_R)
    AT = typeof(at)

    CompositeZiggurat{X,Y,_Ns,M,Mm2,F,ZL,ZM,ZR,AT}(zig_L, zigs_middle, zig_R, at)
end

function RightTailCompositeZiggurat{X,Y,Ns,M}(
    pdf,
    domain;
    ipdfs = nothing,
    area = nothing,
    cdf = nothing,
    ccdf = nothing,
    right_fallback = nothing,
    p = nothing
) where {X,Y,M,Ns<:Tuple{Vararg{Any,M}}}
    _Ns = tupletype_default_numlayers(X, Y, Ns)
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, _Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)

    zig_L = BoundedZiggurat{X,Y,tupletypefirst(_Ns)}(pdf, subdomains[1, :]; ipdf = _ipdfs[1])
    zig_R = UnboundedZiggurat{X,Y,tupletypelast(_Ns)}(
        pdf,
        subdomains[end, :];
        ipdf = _ipdfs[end],
        tailarea = ccdf,
        fallback = right_fallback
    )

    Mm2 = M - 2
    F = typeof(density(zig_L))
    ZL = typeof(zig_L)
    ZM = typeof(zigs_middle)
    ZR = typeof(zig_R)
    AT = typeof(at)

    CompositeZiggurat{X,Y,_Ns,M,Mm2,F,ZL,ZM,ZR,AT}(zig_L, zigs_middle, zig_R, at)
end

function TwoTailCompositeZiggurat{X,Y,Ns,M}(
    pdf,
    domain;
    ipdfs = nothing,
    area = nothing,
    cdf = nothing,
    ccdf = nothing,
    left_fallback = nothing,
    right_fallback = nothing,
    p = nothing
) where {X,Y,M,Ns<:Tuple{Vararg{Any,M}}}
    _Ns = tupletype_default_numlayers(X, Y, Ns)
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, _Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)

    zig_L = UnboundedZiggurat{X,Y,tupletypefirst(_Ns)}(
        pdf,
        subdomains[1, :];
        ipdf = _ipdfs[1],
        tailarea = cdf,
        fallback = left_fallback
    )
    zig_R = UnboundedZiggurat{X,Y,tupletypelast(_Ns)}(
        pdf,
        subdomains[end, :];
        ipdf = _ipdfs[end],
        tailarea = ccdf,
        fallback = right_fallback
    )

    Mm2 = M - 2
    F = typeof(density(zig_L))
    ZL = typeof(zig_L)
    ZM = typeof(zigs_middle)
    ZR = typeof(zig_R)
    AT = typeof(at)

    CompositeZiggurat{X,Y,_Ns,M,Mm2,F,ZL,ZM,ZR,AT}(zig_L, zigs_middle, zig_R, at)
end

function _composite_setup(X, Y, Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)
    dom = regularize(X, domain)

    if M != length(dom) - 1
        error("exactly one number of ziggurat layers (or nothing) must be given for each subdomain, but the number of subdomains does not match the number of layer counts given.")
    end

    if area === nothing
        if cdf !== nothing
            _area = let cdf = cdf
                (a, b) -> abs(cdf(b) - cdf(a))
            end
        elseif ccdf !== nothing
            _area = let ccdf = ccdf
                (a, b) -> abs(ccdf(a) - ccdf(b))
            end
        else
            _tailarea = TailArea{X,Y}(pdf, first(dom))

            _area = let _tailarea=_tailarea
                (a, b) -> abs(_tailarea(b) - _tailarea(a))
            end
        end
    else
        _area = area
    end

    if ipdfs === nothing
        _ipdfs = ntuple(i->nothing, Val(M))
    else
        _ipdfs = ipdfs
    end

    subdomains = get_subdomains(Val(X), Val(M), pdf, dom)

    zigs_middle = middle_zigs(X, Y, Ns, pdf, subdomains, _ipdfs)
    _p = zigprobs(_area, subdomains, p)
    at = AliasTable(_p)

    return zigs_middle, at, _ipdfs, subdomains
end

extracttype(::Type{Type{T}}) where {T} = T
@generated function middle_zigs(X, Y, Ns::Type{<:Tuple}, pdf, subdomains, ipdfs)
    NN = extracttype(Ns).parameters
    idxs = eachindex(NN)

    middles = collect(zip(idxs, NN))[2:(end - 1)]

    Expr(:tuple, (:(BoundedZiggurat{X,Y,$N}(pdf, subdomains[$i, :]; ipdf = ipdfs[$i])) for (i, N) in middles)...)
end

function middle_zigs(X, Y, Ns, pdf, subdomains, ipdfs)
    ntuple(i -> BoundedZiggurat{X,Y,Ns[i + 1]}(pdf, subdomains[i + 1, :]; ipdf = ipdfs[i + 1]), length(Ns)-2)
end

# TODO: add option to autodetect monotonic subdomains.
function composite_ziggurat(pdf, domain; kwargs...)
    domain = regularize(domain)
    Ns = fill(nothing, length(domain) - 1)
    composite_ziggurat(pdf, domain, Ns; kwargs...)
end

function composite_ziggurat(pdf, domain, N::Integer; kwargs...)
    domain = regularize(domain)
    Ns = fill(N, length(domain)-1)
    composite_ziggurat(pdf, domain, Ns; kwargs...)
end

"""
    composite_ziggurat(pdf, domain, Ns; [ipdf, cdf, ccdf, left_fallback, right_fallback, p])

Constructs a sampler for a piecewise-monotonic univariate probability distribution defined by
a probability density function (`pdf`). The domain must be a list of numbers that divides
the `pdf` into monotonic segments. The pdf must be monotonic on each subinterval. It must
not diverge to infinity anywhere on the domain, including at the division points, but may
otherwise be arbitrary - including discontinuous functions. Generate random numbers by
passing the returned ziggurat object to Julia's `rand` or `rand!` functions.

`composite_ziggurat` independently generates monotonic ziggurats for each subdomain. The
leftmost and rightmost ziggurats may be unbounded. Internal subdomains are necessarily bounded.
The [`monotonic_ziggurat`](@ref) function is used internally. The `ipdf` argument is used to
pass a list of inverse pdfs; one for each subdomain. Use the `cdf` and `ccdf` arguments to
override the left and right tail areas respectively. Similarly, use `left_fallback` and
`right_fallback` if necessary.

# Arguments
 - `pdf`: The probability density function of the desired distribution. It must not diverge \
 to infinity anywhere on the `domain`, including the endpoints. It does not need to be \
 normalized, but `ipdf`, `cdf`, and `ccdf` need to have the same normalization as `pdf`.

 - `domain`: a list of numbers that divides the `pdf` into monotonic segments. `domain` may \
 The values will be promoted to the highest float type present, or Float64 if there \
 are no floats. The type that they are promoted to will be the type produced by sampling \
 from the resulting ziggurat. The domain may be unbounded in either or both directions, \
 for example, `(-Inf, 0, Inf)`.

 - `N`: (Optional) The number of layers in each ziggurat. `N` may be a single number that \
 applies to all ziggurats, or a list of numbers one for each subdomain. If `N` is a power \
 of two and the domain is Float64 with N <= 4096, Float32 with N <= 512, or Float16 with \
 N <= 64, then an optimized sampling algorithm is used. Normally, `N` defaults to 256, \
 but for Float16 domains it defaults to 64.

 - `ipdf`: (Optional) A list of functions that invert the pdf on each monotonic subdomain. \
 They are used in the ziggurat construction algorithm, but not during sampling, so they \
 may affect the performance of the initial call to `z = CompositeZiggurat(...)`, but not \
 subsequent calls to `rand(z)`. If `ipdf` is `nothing` or not provided, then all of the \
 ipdfs will be calculated numerically using a root finding algorithm. If any of the entries \
 in the list are `nothing` then those ipdfs will be calculated numerically. See also \
 [`inversepdf`](@ref).

 - `cdf`: (Optional) The cumulative distribution function. It should be normalized the same \
 way the pdf is normalized. Only used if the leftmost subdomain is unbounded below.

 - `ccdf`: (Optional) The complementary cumulative distribution function. It should be \
 normalized the same way the pdf is normalized. Only used if the rightmost subdomain is \
 unbounded above.

 - `left_fallback` and `right_fallback: (Optional) functions that take arguments (rng, a) \
 and produce a random variate in the tail beyond `a` (below `a` for `left_fallback`, and \
 above `a` for `right_fallback`). These are only used on unbounded domains, one or both \
 may be ignored if the domain is bounded on one or both sides. If no fallback is provided, \
 then inverse transform sampling will be used by numerically inverting the `cdf` and/or \
 the `ccdf` (which may themselves be numerical estimates). Overriding the fallbacks can \
 improve the worst-case performance of the sampling algorithm, but since the fallback is \
 only called rarely, it will have a smaller influence on the average performance.

# Examples
```julia-repl
julia> z = CompositeZiggurat(x->exp(-x^2), (-Inf, 0, Inf), 512)
CompositeZiggurat{...}(...)

julia> rand(z, 3)
3-element Vector{Float64}:
  0.756124133407901
 -1.014527969936336
 -1.3362254231230777
```
"""
function composite_ziggurat(
    pdf,
    domain,
    Ns;
    ipdfs = nothing,
    area = nothing,
    cdf = nothing,
    ccdf = nothing,
    left_fallback = nothing,
    right_fallback = nothing,
    p = nothing
)
    dom = regularize(domain)
    X = eltype(dom)
    Y = guess_ytype(pdf, dom)
    Ns = default_numlayers.(X, Y, Ns)
    _Ns = Tuple{Ns...}
    M = length(dom) - 1
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, _Ns, M, pdf, dom; ipdfs, area, cdf, ccdf, p)

    zig_L = monotonic_ziggurat(pdf, subdomains[1, :], Ns[1]; ipdf = _ipdfs[1], cdf = cdf, fallback = left_fallback)
    zig_R =
        monotonic_ziggurat(pdf, subdomains[end, :], Ns[end]; ipdf = _ipdfs[end], ccdf = ccdf, fallback = right_fallback)

    Mm2 = M - 2
    F = typeof(density(zig_L))
    ZL = typeof(zig_L)
    ZM = typeof(zigs_middle)
    ZR = typeof(zig_R)
    AT = typeof(at)

    CompositeZiggurat{X,Y,_Ns,M,Mm2,F,ZL,ZM,ZR,AT}(zig_L, zigs_middle, zig_R, at)
end

"""
    monotonic_segments(f, domain)

Divides the domain into its monotonic subdomains using automatic differentiation and root
finding. The result will have every entry of `domain` as well as the detected critical points
of `f` in sorted order with duplicates dropped.
"""
function monotonic_segments(f, domain)
    domain = regularize(domain)
    df = x -> ForwardDiff.derivative(f, x)
    boundaries = find_zeros(df, domain)
    unique(sort!(vcat(domain.a, boundaries)))
end

"""
    get_subdomains(f, domain)

Convert a list of subdomain boundaries -- `domain` -- into a list of pairs representing subdomains.
Attempts to detect discontinuities in f, and account for them.

# Examples

```julia-repl
julia> Ziggurats.get_subdomains(cos, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, 0.0]
 [0.0, 1.0]

julia> Ziggurats.get_subdomains(x -> x>0 ? 1.0 : 0.0, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, 0.0]
 [5.0e-324, 1.0]

julia> Ziggurats.get_subdomains(x -> x>=0 ? 1.0 : 0.0, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, -5.0e-324]
 [0.0, 1.0]

julia> Ziggurats.get_subdomains(sign, [-1, 0, 1])
2-element Vector{Vector{Float64}}:
 [-1.0, -5.0e-324]
 [5.0e-324, 1.0]
```
"""
function get_subdomains(::Val{X}, ::Val{M}, f, domain::Regularized) where {X,M}
    sd = Array{X}(undef, M, 2)

    sd[1, 1] = domain[1]
    for i in 2:M
        x = domain[i]
        if isapprox(f(prevfloat(x)), f(x); atol = sqrt(eps(x)))
            sd[i - 1, 2] = x
        else
            sd[i - 1, 2] = prevfloat(x)
        end

        if isapprox(f(nextfloat(x)), f(x); atol = sqrt(eps(x)))
            sd[i, 1] = x
        else
            sd[i, 1] = nextfloat(x)
        end
    end
    sd[end, 2] = domain[end]
    sd
end

zigprobs(area, subdomains, ::Nothing) = [area(a, b) for (a, b) in eachrow(subdomains)]
function zigprobs(area, subdomains, p::AbstractArray)
    if length(p) != size(subdomains, 1)
        error("the length of p is not equal to the number of subdomains.")
    end

    _p = similar(p, eltype(subdomains))
    for i in eachindex(p)
        if p[i] === nothing
            a, b = @view(subdomains[i, :])
            _p[i] = area(a, b)
        else
            _p[i] = p[i]
        end
    end
    _p
end

# special case of no just two ziggurats, no 'middle' ziggurats
@inline function Base.rand(
    rng::AbstractRNG,
    s::Random.SamplerTrivial{<:CompositeZiggurat{X,Y,Ns,2,0,F,ZL,Tuple{},ZR}}
) where {X,Y,Ns,F,ZL,ZR}
    cz = s[]

    i = rand(rng, cz.at)
    if i == 1
        return rand(rng, cz.zl)
    else
        return rand(rng, cz.zr)
    end
end

# special case of all 'middle' ziggurats the same
@inline function Base.rand(
    rng::AbstractRNG,
    s::Random.SamplerTrivial{<:CompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR}}
) where {X,Y,Ns,M,Mm2,F,ZL,ZM<:NTuple{Mm2,BoundedZiggurat{X,Y,N,K,NP1,LM,F}} where {N,K,NP1,LM},ZR}
    cz = s[]

    i = rand(rng, cz.at)
    if i == 1
        return rand(rng, cz.zl)
    elseif i == M
        return rand(rng, cz.zr)
    else
        # The type signature of this method guarantees that this indexing operation is type-stable
        return rand(rng, cz.zm[i - 1])
    end
end

# general case
@generated function Base.rand(
    rng::AbstractRNG,
    s::Random.SamplerTrivial{<:CompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR}}
) where {X,Y,Ns,M,Mm2,F,ZL,ZM,ZR}
    # The expression built up here is an if-elseif-else block. Here's an example for M=4:
    #
    # if i == 1
    #     return rand(rng, cz.zl)
    # elseif i == 2
    #     return rand(rng, cz.zm[1]
    # elseif i == 3
    #     return rand(rng, cz.zm[2])
    # else
    #     return rand(rng, cz.zr)
    # end
    #
    # It looks like you could just do `rand(rng, cz.zm[i-1])` and not need a generated function,
    # like the method above, but when the indexing is type-unstable it's very slow.
    # This method is much faster (up to 1000x!). I'm not exactly sure why, but I think it's
    # doing something similar to manual union splitting.

    ex = Expr(:elseif, :(i == $(M-1)), :(return rand(rng, cz.zm[$(M - 1 - 1)])), :(return rand(rng, cz.zr)))
    for j in (M - 2):-1:2
        ex = Expr(:elseif, :(i == $j), :(return rand(rng, cz.zm[$(j - 1)])), ex)
    end
    ex = Expr(:if, :(i == 1), :(return rand(rng, cz.zl)), ex)

    quote
        @inline
        cz = s[]

        i = rand(rng, cz.at)
        $ex
    end
end

@inline function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:CompositeZiggurat{X,Y,Ns,2,0,F,ZL,Tuple{},ZR}}
) where {X<:FloatXX,Y,Ns,F,ZL,ZR}
    cz = s[]
    if length(A) < 7 # TODO: Tune this number
        for i in eachindex(A)
            @inbounds A[i] = rand(rng, s)
        end
    else
        T = corresponding_uint(X)
        GC.@preserve A rand!(rng, Random.UnsafeView{T}(pointer(A), length(A)))

        for i in eachindex(A)
            @inbounds r = reinterpret(T, A[i])

            j = rand(rng, cz.at)
            if j == 1
                w = widths(cz.zl)
                k = layerratios(cz.zl)
                y = heights(cz.zl)
                mb = highside(cz.zl)
                am = lowside(cz.zl)
                pdf = density(cz.zl)
                fb = fallback(cz.zl)
                LM = mask(cz.zl)

                @inbounds A[i] = zigsample2(rng, w, k, y, mb, am, pdf, fb, LM, r)
            else
                w = widths(cz.zr)
                k = layerratios(cz.zr)
                y = heights(cz.zr)
                mb = highside(cz.zr)
                am = lowside(cz.zr)
                pdf = density(cz.zr)
                fb = fallback(cz.zr)
                LM = mask(cz.zr)

                @inbounds A[i] = zigsample2(rng, w, k, y, mb, am, pdf, fb, LM, r)
            end
        end
    end
    A
end

@inline function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:CompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR}}
) where {X<:FloatXX,Y,Ns,M,Mm2,F,ZL,ZM<:NTuple{Mm2,BoundedZiggurat{X,Y,N,K,NP1,LM,F}} where {N,K,NP1,LM},ZR}
    cz = s[]
    if length(A) < 7 # TODO: Tune this number
        for i in eachindex(A)
            @inbounds A[i] = rand(rng, s)
        end
    else
        T = corresponding_uint(X)
        GC.@preserve A rand!(rng, Random.UnsafeView{T}(pointer(A), length(A)))

        for i in eachindex(A)
            @inbounds r = reinterpret(T, A[i])

            j = rand(rng, cz.at)
            if j == 1
                w = widths(cz.zl)
                k = layerratios(cz.zl)
                y = heights(cz.zl)
                mb = highside(cz.zl)
                am = lowside(cz.zl)
                pdf = density(cz.zl)
                fb = fallback(cz.zl)
                LM = mask(cz.zl)

                @inbounds A[i] = zigsample2(rng, w, k, y, mb, am, pdf, fb, LM, r)
            elseif j == M
                w = widths(cz.zr)
                k = layerratios(cz.zr)
                y = heights(cz.zr)
                mb = highside(cz.zr)
                am = lowside(cz.zr)
                pdf = density(cz.zr)
                fb = fallback(cz.zr)
                LM = mask(cz.zr)

                @inbounds A[i] = zigsample2(rng, w, k, y, mb, am, pdf, fb, LM, r)
            else
                z = cz.zm[j - 1]
                w = widths(z)
                k = layerratios(z)
                y = heights(z)
                mb = highside(z)
                am = lowside(z)
                pdf = density(z)
                fb = fallback(z)
                LM = mask(z)

                @inbounds A[i] = zigsample2(rng, w, k, y, mb, am, pdf, fb, LM, r)
            end
        end
    end
    A
end

@generated function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:CompositeZiggurat{X,Y,Ns,M,Mm2,F,ZL,ZM,ZR}}
) where {X<:FloatXX,Y,Ns,M,Mm2,F,ZL,ZM,ZR}
    sample_ex = z_ex -> quote
        w = widths($z_ex)
        k = layerratios($z_ex)
        y = heights($z_ex)
        mb = highside($z_ex)
        am = lowside($z_ex)
        pdf = density($z_ex)
        fb = fallback($z_ex)
        LM = mask($z_ex)

        @inbounds A[i] = zigsample2(rng, w, k, y, mb, am, pdf, fb, LM, r)
    end

    ex = Expr(:elseif, :(i == $(M-1)), sample_ex(:(cz.zm[$(M - 1 - 1)])), sample_ex(:(cz.zr)))
    for j in (M - 2):-1:2
        ex = Expr(:elseif, :(i == $j), sample_ex(:(cz.zm[$(j - 1)])), ex)
    end
    ex = Expr(:if, :(i == 1), sample_ex(:(cz.zl)), ex)

    quote
        cz = s[]
        if length(A) < 7 # TODO: Tune this number
            for i in eachindex(A)
                @inbounds A[i] = rand(rng, s)
            end
        else
            T = corresponding_uint(X)
            GC.@preserve A rand!(rng, Random.UnsafeView{T}(pointer(A), length(A)))

            for i in eachindex(A)
                @inbounds r = reinterpret(T, A[i])

                j = rand(rng, cz.at)

                $ex
            end
        end
        A
    end
end
