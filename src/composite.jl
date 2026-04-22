struct CompositeZiggurat{
    X,
    Y,
    Ns,
    M,
    ZL<:MonotonicZiggurat{X,Y},
    ZR<:MonotonicZiggurat{X,Y},
    Z<:NTuple{M,MonotonicZiggurat{X,Y}},
    AT
} <: Ziggurat{X,Y}
    zigs::Z
    at::AT

    function CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT}(zigs, at) where {X,Y,Ns<:Tuple,M,ZL,ZR,Z,AT}
        if Z.parameters[1] != ZL || Z.parameters[end] != ZR
            error("the first and last entries of Z must match ZL and ZR respectively.")
        end

        for (z, n) in zip(Z.parameters, Ns.parameters)
            if numlayers(z) != n
                error("The number of layers in each entry of Z must equal the cooresponding entry of Ns")
            end
        end

        for z in Z.parameters[2:(end - 1)]
            if !(z <: BoundedZiggurat)
                error("each entry in Z except the first and last must be BoundedZiggurat.")
            end
        end

        new{X,Y,Ns,M,ZL,ZR,Z,AT}(zigs, at)
    end
end

const BoundedCompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT} =
    CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT} where {X,Y,Ns,M,ZL<:BoundedZiggurat,ZR<:BoundedZiggurat,Z,AT}
const LeftTailCompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT} =
    CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT} where {X,Y,Ns,M,ZL<:UnboundedZiggurat,ZR<:BoundedZiggurat,Z,AT}
const RightTailCompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT} =
    CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT} where {X,Y,Ns,M,ZL<:BoundedZiggurat,ZR<:UnboundedZiggurat,Z,AT}
const TwoTailCompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT} =
    CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT} where {X,Y,Ns,M,ZL<:UnboundedZiggurat,ZR<:UnboundedZiggurat,Z,AT}

tupletypelength(::Type{T}) where {T<:Tuple} = fieldcount(T)

tupletypefirst(::Type{<:Tuple{First,Vararg}}) where {First} = First
tupletypelast(::Type{T}) where {T<:Tuple{Any,Vararg}} = fieldtype(T, fieldcount(T))

asatupletype(t::Type{<:Tuple}) = t
asatupletype(t) = Tuple{t...}

handle_N_as_value(N::Integer, M) = NTuple{M,N}
handle_N_as_value(N::Tuple, M) = asatupletype(N)

for CZ in
    (:BoundedCompositeZiggurat, :LeftTailCompositeZiggurat, :RightTailCompositeZiggurat, :TwoTailCompositeZiggurat)
    @eval begin
        function $(CZ){X,Y,N,M}(pdf, domain; kwargs...) where {X,Y,N,M}
            _N = handle_N_as_value(N, M)
            $(CZ){X,Y,_N,M}(pdf, domain; kwargs...)
        end

        function $(CZ){X,Y,Ns}(pdf, domain; kwargs...) where {X,Y,Ns}
            _Ns = asatupletype(Ns)
            M = tupletypelength(_Ns)
            $(CZ){X,Y,_Ns,M}(pdf, domain; kwargs...)
        end

        $(CZ){X,Y,Ns,M}(pdf, domain; kwargs...) where {X,Y,Ns<:Tuple,M} =
            error("the Tuple of layer counts, Ns, must have length equal to M.")
    end
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
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)

    zig_L = BoundedZiggurat{X,Y,tupletypefirst(Ns)}(pdf, subdomains[1, :]; ipdf = _ipdfs[1])
    zig_R = BoundedZiggurat{X,Y,tupletypelast(Ns)}(pdf, subdomains[end, :]; ipdf = _ipdfs[end])

    zigs = (zig_L, zigs_middle..., zig_R)

    ZL = typeof(zig_L)
    ZR = typeof(zig_R)
    Z = typeof(zigs)
    AT = typeof(at)

    CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT}(zigs, at)
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
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)

    zig_L = UnboundedZiggurat{X,Y,tupletypefirst(Ns)}(
        pdf,
        subdomains[1, :];
        ipdf = _ipdfs[1],
        tailarea = cdf,
        fallback = left_fallback
    )
    zig_R = BoundedZiggurat{X,Y,tupletypelast(Ns)}(pdf, subdomains[end, :]; ipdf = _ipdfs[end])

    zigs = (zig_L, zigs_middle..., zig_R)

    ZL = typeof(zig_L)
    ZR = typeof(zig_R)
    Z = typeof(zigs)
    AT = typeof(at)

    CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT}(zigs, at)
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
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)

    zig_L = BoundedZiggurat{X,Y,tupletypefirst(Ns)}(pdf, subdomains[1, :]; ipdf = _ipdfs[1])
    zig_R = UnboundedZiggurat{X,Y,tupletypelast(Ns)}(
        pdf,
        subdomains[end, :];
        ipdf = _ipdfs[end],
        tailarea = ccdf,
        fallback = right_fallback
    )

    zigs = (zig_L, zigs_middle..., zig_R)

    ZL = typeof(zig_L)
    ZR = typeof(zig_R)
    Z = typeof(zigs)
    AT = typeof(at)

    CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT}(zigs, at)
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
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)

    zig_L = UnboundedZiggurat{X,Y,tupletypefirst(Ns)}(
        pdf,
        subdomains[1, :];
        ipdf = _ipdfs[1],
        tailarea = cdf,
        fallback = left_fallback
    )
    zig_R = UnboundedZiggurat{X,Y,tupletypelast(Ns)}(
        pdf,
        subdomains[end, :];
        ipdf = _ipdfs[end],
        tailarea = ccdf,
        fallback = right_fallback
    )

    zigs = (zig_L, zigs_middle..., zig_R)

    ZL = typeof(zig_L)
    ZR = typeof(zig_R)
    Z = typeof(zigs)
    AT = typeof(at)

    CompositeZiggurat{X,Y,Ns,M,ZL,ZR,Z,AT}(zigs, at)
end

function _composite_setup(X, Y, Ns, M, pdf, domain; ipdfs, area, cdf, ccdf, p)
    dom = regularize(domain)

    if M != length(dom) - 1
        error("the number of subdomains must match the number of layernums.")
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
            segbuf = alloc_segbuf(X, Y, Y)
            _tailarea = TailArea(pdf, first(dom), segbuf)

            _area = let _tailarea=_tailarea
                (a, b) -> abs(_tailarea(b) - _tailarea(a))
            end
        end
    else
        _area = area
    end

    if ipdfs === nothing
        _ipdfs = ntuple(i->nothing, M)
    else
        _ipdfs = ipdfs
    end

    subdomains = get_subdomains(pdf, dom)

    zigs_middle = middle_zigs(X, Y, Ns, pdf, subdomains, _ipdfs)
    _p = zigprobs(_area, subdomains, p)
    at = AliasTable(_p)

    return zigs_middle, at, _ipdfs, subdomains
end

extracttype(::Type{Type{T}}) where {T} = T
@generated function middle_zigs(X, Y, Ns::Type{<:Tuple}, pdf, subdomains, ipdfs)
    # TODO: Check if this @generated method is needed for type stability. The generic method
    # (below) may be sufficient.
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
    N = default_numlayers(nothing, eltype(domain))
    composite_ziggurat(pdf, domain, N; kwargs...)
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
    M = length(dom) - 1
    zigs_middle, at, _ipdfs, subdomains = _composite_setup(X, Y, Ns, M, pdf, dom; ipdfs, area, cdf, ccdf, p)

    zig_L = monotonic_ziggurat(pdf, subdomains[1, :], Ns[1]; ipdf = _ipdfs[1], cdf = cdf, fallback = left_fallback)
    zig_R =
        monotonic_ziggurat(pdf, subdomains[end, :], Ns[end]; ipdf = _ipdfs[end], ccdf = ccdf, fallback = right_fallback)

    zigs = (zig_L, zigs_middle..., zig_R)

    ZL = typeof(zig_L)
    ZR = typeof(zig_R)
    Z = typeof(zigs)
    AT = typeof(at)

    CompositeZiggurat{X,Y,asatupletype(Ns),M,ZL,ZR,Z,AT}(zigs, at)
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
get_subdomains(f, domain) = _get_subdomains(f, regularize(domain))
get_subdomains(f, domain::Regularized) = _get_subdomains(f, domain)
function _get_subdomains(f, domain)
    sd = Array{eltype(domain)}(undef, length(domain)-1, 2)

    sd[1, 1] = domain[1]
    for (i, x) in enumerate(domain[2:(end - 1)])
        if isapprox(f(prevfloat(x)), f(x); atol = sqrt(eps(x)))
            sd[i, 2] = x
        else
            sd[i, 2] = prevfloat(x)
        end

        if isapprox(f(nextfloat(x)), f(x); atol = sqrt(eps(x)))
            sd[i + 1, 1] = x
        else
            sd[i + 1, 1] = nextfloat(x)
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

@inline function Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:CompositeZiggurat})
    cz = s[]
    zigs = cz.zigs
    i = rand(rng, cz.at)
    @inbounds z = zigs[i]
    rand(rng, z)
end

@inline function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:CompositeZiggurat{X}}
) where {X<:FloatXX}
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
            @inbounds z = cz.zigs[j]

            w = widths(z)
            k = layerratios(z)
            y = heights(z)
            mb = highside(z)
            am = lowside(z)
            pdf = density(z)
            fb = fallback(z)
            LM = mask(z)

            @inbounds A[i] = zigsample(rng, w, k, y, mb, am, pdf, fb, LM, r)
        end
    end
    A
end
