# Type parameters: domain type, range type, Layer Mask (unsigned integer), Rearrange layer bits (Bool).
abstract type MonotonicZiggurat{X,Y,LM} <: Ziggurat{X,Y} end

struct BoundedZiggurat{X,Y,LM,K,N,NP1,F} <: MonotonicZiggurat{X,Y,LM}
    w::SVector{NP1,X}
    k::SVector{N,K}
    y::SVector{NP1,Y}
    modalboundary::X
    pdf::F
    function BoundedZiggurat(w, k, y, mb, pdf)
        LM = layermask(eltype(w), length(k))
        new{eltype(w),eltype(y),LM,eltype(k),length(k),length(w),typeof(pdf)}(w, k, y, mb, pdf)
    end
end

struct UnboundedZiggurat{X,Y,LM,K,N,NP1,F,FB} <: MonotonicZiggurat{X,Y,LM}
    w::SVector{NP1,X}
    k::SVector{N,K}
    y::SVector{NP1,Y}
    modalboundary::X
    pdf::F
    fallback::FB
    function UnboundedZiggurat(w, k, y, mb, pdf, fb)
        LM = layermask(eltype(w), length(k))
        new{eltype(w),eltype(y),LM,eltype(k),length(k),length(w),typeof(pdf),typeof(fb)}(w, k, y, mb, pdf, fb)
    end
end

widths(z::Ziggurat) = z.w
numlayers(z::Ziggurat) = length(widths(z)) - 1
layerratios(z::Ziggurat) = z.k
heights(z::Ziggurat) = z.y
highside(z::Ziggurat) = z.modalboundary
density(z::Ziggurat) = z.pdf
fallback(::BoundedZiggurat) = nothing
fallback(z::UnboundedZiggurat) = z.fallback
xarray(z::Ziggurat) = significand_bitmask(eltype(z)) .* widths(z) .+ highside(z)

corresponding_uint(::Type{Float64}) = UInt64
corresponding_uint(::Type{Float32}) = UInt32
corresponding_uint(::Type{Float16}) = UInt16

const FloatXX = Union{Float64,Float32,Float16}

significand_bits(::Type{Float64}) = 52
significand_bits(::Type{Float32}) = 23
significand_bits(::Type{Float16}) = 10

# shiftbits is the number of bits besides the significand (i.e. exponent bits plus one)
shiftbits(::Type{Float64}) = 12
shiftbits(::Type{Float32}) = 9
shiftbits(::Type{Float16}) = 6

function layer_bits(T, LM, r)
    UT = corresponding_uint(T)
    u_mask = significand_bitmask(T) << shiftbits(T)
    overlap_mask = LM & u_mask
    if overlap_mask == 0
        return r & LM
    end
    other_mask = UT(2)^shiftbits(T) - UT(1) # first `shiftbits(T)` bits set
    rearrange_shift = ilog2(overlap_mask >>> shiftbits(T)) + 1

    overlap_bits = (r & overlap_mask) >>> shiftbits(T)
    other_bits = (r & other_mask) << rearrange_shift
    overlap_bits | other_bits
end

function layer_bits_signed(T, LM, r)
    UT = corresponding_uint(T)
    u_mask = significand_bitmask(T) << shiftbits(T)
    overlap_mask = LM & u_mask
    if overlap_mask == 0
        return (r & LM) >>> 1
    end
    other_mask = (UT(2)^shiftbits(T) - UT(1)) - UT(1) # first `shiftbits(T)` bits set except the first one
    rearrange_shift = ilog2(overlap_mask >>> shiftbits(T))

    overlap_bits = (r & overlap_mask) >>> shiftbits(T)
    other_bits = (r & other_mask) << rearrange_shift
    overlap_bits | other_bits
end

function layermask(T, N)
    if !isinteger(N) || N < 1
        throw(ArgumentError("N must be an integer greater than zero."))
    end
    _layermask(T, N)
end
_layermask(::Type, N) = nothing, false
function _layermask(T::Type{<:FloatXX}, N)
    UT = corresponding_uint(T)
    if !ispow2(N) || N > big(typemax(UT)) + 1
        return nothing, false
    end
    UT(N - 1)
end

"""
    layermask with one bit reserved for the sign
"""
function layermask_signed(T, N)
    if !isinteger(N) || N < 1
        throw(ArgumentError("N must be an integer greater than zero."))
    end
    _layermask_signed(T, N)
end
_layermask_signed(::Type, N) = nothing
function _layermask_signed(T::Type{<:FloatXX}, N)
    UT = corresponding_uint(T)
    if !ispow2(N) || N > UInt64(2)^(8sizeof(T) - 1)
        return nothing
    end
    UT(N - 1) << 1
end

"""
Provide error messages for unexpected or invalid results including non-montonicity.
"""
struct PDFWrap{F,X,Y}
    f::F
    mb::X
    am::X
    fmb::Y
    fam::Y
    function PDFWrap(f, mb, am)
        mb, am = promote(mb, am)
        fmb = f(mb)
        fam = f(am)
        fmb, fam = promote(fmb, fam)
        if isnan(fmb) || isnan(fam)
            error("pdf is NaN on the boundary. Check the definition of your pdf. It must \
            return non-negative numbers everywhere on its domain, including the end points \
            of the domain.")
        end
        if fmb < fam
            error("modalboundary and argminboundary are misidentified.")
        end
        if fam < zero(fam)
            error("pdf($am) is negative. Check the definition of your pdf. It must return \
            non-negative numbers everywhere on its domain, including the end points of \
            the domain.")
        end
        if fam == fmb == zero(fam)
            error("pdf is zero on both endpoints of the domain, $(minmax(am,mb)). Ensure \
            that the pdf is monotonic.")
        end
        new{typeof(f),typeof(mb),typeof(fmb)}(f, mb, am, fmb, fam)
    end
end

function (pdf::PDFWrap)(x)
    if !between(pdf.mb, pdf.am, x)
        error("Unexpected Error: attempted to evaluate pdf at x = $x outside the domain.
        Please report this error.")
    end

    result = pdf.f(x)

    if isnan(result)
        error("pdf($x) is NaN. Check the definition of your pdf. It must \
        return positive numbers everywhere on its domain, including the end points of \
        the domain.")
    end

    if isinf(result)
        error("pdf($x) is infinite. The Marsaglia & Tsang Ziggurat \
        algorithm requires a finite pdf.")
    end

    δ = √eps(pdf.fmb)
    if !between(pdf.fmb + δ, pdf.fam - δ, result)
        d = minmax(pdf.mb, pdf.am)
        error("pdf is not monotonic on the domain = ($(d[1]), $(d[2])), pdf($x) = $result.")
    end

    return float(result)
end

"""
Force the ipdf to respect the domain and codomain.
"""
struct IPDFWrap{F,X,Y}
    ipdf::F
    mb::X
    am::X
    fmb::Y
    fam::Y
    function IPDFWrap(ipdf, mb, am, fmb, fam)
        mb, am = promote(mb, am)
        fmb, fam = promote(fmb, fam)
        new{typeof(ipdf),typeof(mb),typeof(fmb)}(ipdf, mb, am, fmb, fam)
    end
end

# TODO: make inverse and IPDFWrap give the same error messages.
function (ipdf::IPDFWrap)(y)
    if y < zero(y)
        error("Unexpected error: ipdf was about to be evaluated at y = $y, but the argument \
        to ipdf should always be positive.")
    elseif y <= ipdf.fam
        # handle the discontinuity at the boundary
        return ipdf.am
    elseif y == ipdf.fmb
        # Ideally, ipdf(f(mb)) = mb would be true by definition, but floating point
        # inaccuracies can cause problems. Making this evaluation exact may help avoid
        # domain errors.
        return ipdf.mb
    elseif y > ipdf.fmb
        error("ipdf evaluated at y = $y which is greater than the pdf at either boundary. \
        Make sure your pdf is monotonic and the inverse pdf is correct.")
    else
        result = ipdf.ipdf(y)
        lb, ub = minmax(ipdf.mb, ipdf.am)
        return clamp(result, lb, ub)
    end
end

"""
    NoWrap(pdf)
    NoWrap(ipdf)

Forces the ziggurat construction to use unwrapped functions internally.

The ziggurat constructors wrap the pdf and ipdf functions to guarentee that they respect
certain assumptions. E.g., the pdf is non-negative and finite, and the ipdf respects the
domain given to the constructor. If you pass the functions with NoWrap it forces the ziggurat
constructor to use the unwrapped function throughout.

While this function is part of the public interface, it does interact with non-public internals.
It is provided for the unlikely event that the wrapper is making an error.
"""
struct NoWrap{F}
    f::F
end

PDFWrap(f::NoWrap, mb, am, fmb, fam) = f.f
IPDFWrap(f::NoWrap, mb, am, fmb, fam) = f.f

# TODO: Need tests for unwrapped functions.
"""
    monotonic_ziggurat(pdf, domain, [N]; [ipdf, tailarea, fallback, ...])

Constructs a high-performance sampler for a univariate probability distribution defined by a
probability density function (`pdf`). The pdf must be monotonic on the domain and must not
diverge to infinity anywhere on the domain, including at the endpoints, but may otherwise be
arbitrary - including discontinuous functions. Generate random numbers by passing the returned
ziggurat object to Julia's `rand` or `rand!` functions.

`monotonic_ziggurat` selects between constructing a [`BoundedZiggurat`](@ref) or an
[`UnboundedZiggurat`](@ref) based on whether the domain is finite. Both types of ziggurats
use an inverse pdf in their construction. Normally it is computed using a root finding method,
but it may be overridden using the `ipdf` argument. `UnboundedZiggurat`s also use a `tailarea`
function during construction and a `fallback` during sampling. Normally these additional
functions are computed numerically, but they can be provided explicitly as keyword arguments
if necessary.

# Arguments
 - `pdf`: The probability density function of the desired distribution. It must be monotonic \
 and must not diverge to infinity anywhere on the `domain`, including the endpoints. It \
 does not need to be normalized, but `ipdf` and `tailarea` need to have the same \
 normalization as `pdf`.

 - `domain`: The domain of the pdf. `domain` may be any collection of numbers, but only its \
 extrema will be used as the boundaries of the domain. The values will be promoted to the \
 highest float type present, or Float64 if there are no floats. The type that they are \
 promoted to will be the type produced by sampling from the resulting ziggurat. The domain \
 may be unbounded, for example, `(0, Inf)`.

 - `N`: (Optional) The number of layers in the ziggurat. If `N` is a power of two and the \
 domain is Float64 with N <= 4096, Float32 with N <= 512, or Float16 with N <= 64, then \
 an optimized sampling algorithm is used. Normally, `N` defaults to 256, but for Float16 \
 domains it defaults to 64.

 - `ipdf`: (Optional) This function is the inverse of the given `pdf` argument. It's used in \
 the ziggurat construction algorithm, but not during sampling, so it may affect the \
 performance of the initial call to `z = monotonic_ziggurat(...)`, but not subsequent \
 calls to `rand(z)`. If no ipdf is provided then the inverse will be calculated numerically \
 using a root finding algorithm. See also [`inversepdf`](@ref).

 - `tailarea`: (Optional) A function that computes the area of the tail of the `pdf` when the \
 domain is unbounded. If the domain has a finite length, this argument is ignored. The \
 `tailarea` argument takes precedence over the `cdf` and `ccdf` arguments. If neither \
 `tailarea`, nor an appropriate (c)cdf argument is provided then a numerical integral \
 will be computed. Overriding `tailarea` may have an indirect effect on sampling performance \
 since the default fallback algorithm uses `tailarea`.

 - `cdf`: (Optional) The cumulative distribution function. It should be normalized the same \
 way the pdf is normalized. Superceeded by `tailarea`, and ignored when the domain is \
 bounded.

 - `ccdf`: (Optional) The complementary cumulative distribution function. It should be \
 normalized the same way the pdf is normalized. Superceeded by `tailarea`, and ignored \
 when the domain is bounded.

 - `fallback`: (Optional) A function that takes arguments (rng, a) and produces a random \
 variate in the tail beyond `a` (below `a` for increasing distributions, or above `a` \
 for decreasing distributions). This is only used on unbounded domains, it is ignored on \
 bounded domains. If no fallback is provided, then inverse transform sampling will be used \
 by numerically inverting the `tailarea` (which may itself may be a numerical estimate). \
 Overriding the fallback can greatly improve the worst-case performance of the sampling \
 algorithm, but since the fallback is only called rarely, it will have a smaller influence \
 on the average performance.

# Examples
```julia-repl
julia> z = monotonic_ziggurat(x->exp(-x), (0,Inf))
UnboundedZiggurat{...}(...)

julia> rand(z)
0.3594408084987401

julia> z = monotonic_ziggurat(sin, (0,π/2), 512)
BoundedZiggurat{...}(...)

julia> rand(z, 100)
100-element Vector{Float64}
[...]

julia> z = monotonic_ziggurat(x->exp(-x^5), (0.0f0, Inf32), 512)
UnboundedZiggurat{...}(...)

julia> rand(z)
0.29527622f0
```
"""
function monotonic_ziggurat(
    pdf,
    domain,
    N = nothing;
    ipdf = inversepdf(pdf, domain),
    tailarea = nothing,
    cdf = nothing,
    ccdf = nothing,
    fallback = nothing
)
    domain = regularize(domain)

    if N === nothing
        N = eltype(domain) == Float16 ? 64 : 256
    end

    a, b = extrema(domain)
    if isinf(a) || isinf(b)
        tailarea = _choose_tailarea_func(pdf, domain, tailarea, cdf, ccdf)
        UnboundedZiggurat(pdf, domain, N; ipdf, tailarea, fallback)
    else
        BoundedZiggurat(pdf, domain, N; ipdf)
    end
end

"""
    BoundedZiggurat(pdf, domain, N; [ipdf])

Constructs a high-performance sampler for a univariate probability distribution defined by a
probability density function (`pdf`) with bounded support. The pdf must be monotonic on the 
domain and must not diverge to infinity anywhere on the domain, including at the endpoints,
but may otherwise be arbitrary - including discontinuous functions. An inverse function to the
pdf, `ipdf`, is needed to construct the sampler. By default, the inverse is computed
numerically, but it can also be provided explicity if necessary. 
"""
function BoundedZiggurat(pdf, domain, N; ipdf = nothing)
    (; x, y, modalboundary) = _bounded_zig_data(pdf, domain, N, ipdf)

    w, k = _optimized_tables(x, modalboundary)

    # final ziggurat uses the unwrapped function so that there is no effect on
    # sampling performance
    BoundedZiggurat(w, k, y, modalboundary, pdf)
end

"""
    UnboundedZiggurat(pdf, domain, N; [ipdf, tailarea, fallback])

Constructs a high-performance sampler for a univariate probability distribution defined by a
probability density function (`pdf`) with unbounded support. The pdf must be monotonic on the
domain and must not diverge to infinity anywhere on the domain, including at the endpoints,
but may otherwise be arbitrary - including discontinuous functions. An inverse pdf and a
tailarea function are used in the construction of the ziggurat, and a `fallback` is used
during sampling. Normally these additional functions are computed numerically, but they can
be provided explicity as keyword arguments if necessary.
"""
function UnboundedZiggurat(pdf, domain, N; ipdf = nothing, tailarea = nothing, fallback = nothing)
    domain = extrema(regularize(domain))
    _check_arguments(N, domain)
    modalboundary, argminboundary = _identify_mode(pdf, domain)

    (; x, y, tailarea) = _unbounded_zig_data(pdf, domain, N, ipdf, tailarea, modalboundary, argminboundary)

    w, k = _optimized_tables(x, modalboundary)

    if fallback === nothing
        # TODO: fallback should come from a 'tool' with this as default but also allows customization.
        # e.g. pass arguments through inverse to find_zero. Think about reducing layers of indirection.
        ta = tailarea(x[2])
        td = minmax(argminboundary, x[2])
        inverse_tailprob = let tailarea = tailarea, ta = ta, td = td
            inversepdf(xx -> tailarea(xx) / ta, td)
        end

        fallback = (rng, _) -> inverse_tailprob(rand(rng, typeof(modalboundary)))
    end

    # final ziggurat uses the unwrapped function so that there is no effect on
    # sampling performance
    UnboundedZiggurat(w, k, y, modalboundary, pdf, fallback)
end

function _bounded_zig_data(pdf, domain, N, ipdf)
    domain = extrema(regularize(domain))

    _check_arguments(N, domain)
    modalboundary, argminboundary = _identify_mode(pdf, domain)

    wpdf = PDFWrap(pdf, modalboundary, argminboundary)

    if ipdf === nothing
        ipdf = inversepdf(wpdf, domain)
    end

    wipdf = IPDFWrap(ipdf, modalboundary, argminboundary, wpdf(modalboundary), wpdf(argminboundary))

    if wpdf(modalboundary) == 0
        error("expected the pdf to be non-zero on at least one boundary.")
    end

    if isinf(domain[1]) || isinf(domain[2])
        error("expected a bounded domain, got domain=$domain.")
    end

    # Build ziggurats using wrapped functions
    x, y = search(N, modalboundary, argminboundary, wpdf, wipdf)

    (; x, y, modalboundary)
end

function _unbounded_zig_data(pdf, domain, N, ipdf, tailarea, modalboundary, argminboundary)
    wpdf = PDFWrap(pdf, modalboundary, argminboundary)

    if ipdf === nothing
        ipdf = inversepdf(wpdf, domain)
    end

    wipdf = IPDFWrap(ipdf, modalboundary, argminboundary, pdf(modalboundary), pdf(argminboundary))

    if wpdf(modalboundary) == 0
        error("expected the pdf to be non-zero on at least one boundary.")
    end

    if !isinf(domain[1]) && !isinf(domain[2])
        error("expected an unbounded domain, got domain=$domain.")
    end

    if tailarea === nothing
        modepdf = wpdf(modalboundary)
        domain_type = typeof(modalboundary)
        range_type = typeof(modepdf)
        error_type = typeof(norm(modepdf))
        segbuf = alloc_segbuf(domain_type, range_type, error_type)

        tailarea = TailArea(wpdf, argminboundary, segbuf)
    end

    # Build ziggurats using wrapped functions
    x, y = search(N, modalboundary, argminboundary, wpdf, wipdf, tailarea)

    (; x, y, tailarea)
end

function _optimized_tables(x, modalboundary)
    w = (x .- modalboundary) ./ significand_bitmask(eltype(x))
    k = [fixedbit_fraction((x[i + 1] - modalboundary)/(x[i] - modalboundary)) for i in 1:(length(x) - 1)]

    w, k
end

struct TailArea{F,X,S}
    f::F
    r::X
    segbuf::S
end

# TODO: a way to pass through arguments to quadgk.
function (ta::TailArea)(x)
    if x == ta.r
        # quadgk returns nan when the domain is empty like that
        # this is a workaround
        zero(ta.r)
    else
        # TODO: Add error tolerance and check the returned error estimate.
        abs(quadgk(ta.f, x, ta.r; ta.segbuf)[1])
    end
end

function _choose_tailarea_func(pdf, domain, tailarea, cdf, ccdf)
    if tailarea !== nothing || (cdf === nothing && ccdf === nothing)
        return tailarea
    end

    # below here, tailarea is nothing, and cdf and ccdf are not both nothing.
    a, b = extrema(domain)
    fa, fb = pdf(a), pdf(b)
    if fa < fb # pdf is increasing (constant functions are decreasing)
        if cdf !== nothing
            return cdf
        elseif ccdf !== nothing
            throw(ArgumentError("a ccdf is provided for an increasing pdf, pass cdf or \
            tailarea instead or in addition."))
        else
            error("Unreachable error: it should be impossible to throw this error.")
        end
    else
        if ccdf !== nothing
            return ccdf
        elseif cdf !== nothing
            throw(ArgumentError("a cdf is provided for an increasing pdf, pass ccdf or \
            tailarea instead or in addition."))
        else
            error("Unreachable error: it should be impossible to throw this error.")
        end
    end

    error("Unreachable error: it should be impossible to throw this error.")
end

function _check_arguments(N, domain)
    if N < 1
        throw(DomainError(N, "N must be a positive integer, got N=$N."))
    end
    if isempty(domain)
        error("empty domains are not allowed, got domain=$domain.")
    end

    a, b = extrema(domain)
    if a == b
        error("empty domains are not allowed, got domain=$domain.")
    end
    if isinf(a) && isinf(b)
        error("a domain of (-Inf, Inf) is impossible for a monotonic distribution.")
    end

    return nothing
end

function _identify_mode(pdf, domain)
    # Return the modalboundary (mb) and argminboundary (am).
    # Assume that the domain is well formed and appropriate for a monotonic
    # distribution. A constant function is treated as decreasing.
    a, b = extrema(domain)
    if isinf(a)
        mb = b
        am = a
    elseif isinf(b)
        mb = a
        am = b
    else
        fa, fb = pdf(a), pdf(b)
        if isnan(fa) || isnan(fb)
            error("pdf(x) is NaN on the boundary. Check the definition of your pdf. It must \
            return positive numbers everywhere on its domain, including the end points of \
            the domain.")
        end
        if isinf(fa) || isinf(fb)
            error("pdf(x) is infinite on the boundary. The Marsaglia & Tsang Ziggurat \
            algorithm requires a finite pdf.")
        end
        if fa >= fb
            mb = a
            am = b
        else
            mb = b
            am = a
        end
    end

    mb, am
end

# Bounded
function search(N, modalboundary, argminboundary, pdf, ipdf)
    x, y, y_domain, modalpdf = _search_setup(N, modalboundary, pdf)
    buildargs = (modalboundary, argminboundary, ipdf, modalpdf)
    p = (; x, y, modalpdf, buildargs)
    _search(y_domain, p)
end

# Unbounded
function search(N, modalboundary, argminboundary, pdf, ipdf, tailarea)
    x, y, y_domain, modalpdf = _search_setup(N, modalboundary, pdf)
    buildargs = (modalboundary, argminboundary, ipdf, modalpdf, tailarea)
    p = (; x, y, modalpdf, buildargs)
    _search(y_domain, p)
end

function _search_setup(N, modalboundary, pdf)
    modalpdf = pdf(modalboundary)
    x = Vector{typeof(float(modalboundary))}(undef, N + 1)
    y = Vector{typeof(modalpdf)}(undef, N + 1)

    y_domain = (nextfloat(zero(modalpdf)), modalpdf)

    x, y, y_domain, modalpdf
end

function _search(y_domain, p)
    ystar = find_zero(ziggurat_residual, y_domain, Roots.Bisection(), p)

    x, y = build!(p.x, p.y, ystar, p.buildargs...)

    if isinf(y[end])
        error("ziggurat failed to converge.")
    end

    if y[end] < p.modalpdf
        # Assuming the final interval is [ystar, nextfloat(ystar)], which is true for the Bisection method over
        # Float64, Float32, or Float16's.
        x = copy(x)
        y = copy(y)
        xx, yy = build!(p.x, p.y, nextfloat(ystar), p.buildargs...)

        if isinf(yy[end])
            @warn "the ziggurat doesn't fully cover the pdf, but a larger ziggurat can't be constructed. \
            This is sometimes caused by using too many layers for the level of floating point precision. \
            The ziggurat may be accurate enough to be usable if the height of the ziggurat, \
            `Ziggurats.heights(zig)[end])`, is close to the height of the pdf, `Ziggurats.density(zig, Ziggurats.highside(zig))`."

            return x, y
        end

        if yy[end] < p.modalpdf
            error("the ziggurat search didn't converge as expected.")
        end

        return xx, yy
    end

    return x, y
end

function ziggurat_residual(y2, p)
    _x, _y = build!(p.x, p.y, y2, p.buildargs...)
    _y[end] - p.modalpdf
end

# The caller of build!() is responsible for ensuring consistancy of the build
# arguments. I.e.,
# 1) length(x) == length(y) >= 1
# 2) zero(y2) < y2 <= pdf(boundary)
# 3) boundarypdf = pdf(boundary)
# 4) the pdf is monotonic
# 5) the ipdf is the gerneralized inverse of pdf (i.e. ipdf(y) is the largest x
#       such that pdf(x) >= y for decreasing pdfs, and ipdf(y) is the smallest x
#       such that pdf(x) >= y for increasing pdfs)

# Bounded support
function build!(x, y, y2, modalboundary, argminboundary, ipdf, modalpdf)
    initialize!(y, y2)

    A = layerarea(y[2], modalboundary, argminboundary)
    if A == 0
        # y2 is too small
        y[end] = zero(eltype(y))
        return x, y
    end
    x[1] = argminboundary

    finalize!(x, y, modalboundary, A, ipdf, modalpdf)
end

# Unbounded support
function build!(x, y, y2, modalboundary, argminboundary, ipdf, modalpdf, tailarea)
    initialize!(y, y2)

    x2 = ipdf(y2)
    if isinf(x2)
        # y2 is too small
        y[end] = zero(eltype(y))
        return x, y
    end
    A = layerarea(y[2], x2, modalboundary, tailarea)

    slopesign = sign(modalboundary - argminboundary)
    x[1] = modalboundary - slopesign * A / y[2]

    finalize!(x, y, modalboundary, A, ipdf, modalpdf)
end

function initialize!(y, y2)
    y[1] = zero(eltype(y))
    y[2] = y2
end

# Bounded support
layerarea(y2, modalboundary, argminboundary) = abs(argminboundary - modalboundary) * y2

# Unbounded support
function layerarea(y2, x2, modalboundary, tailarea)
    if y2 == zero(y2)
        return zero(x2 * y2)
    end
    abs(x2 - modalboundary) * y2 + tailarea(x2)
end

function finalize!(x, y, modalboundary, A, ipdf, modalpdf)
    for i in eachindex(x)[(begin + 1):(end - 1)]
        if y[i] >= modalpdf
            # failed to build ziggurat, y2 is too large.
            y[end] = typemax(eltype(y))
            return x, y
        end
        x[i] = ipdf(y[i])
        y[i + 1] = A / abs(x[i] - modalboundary) + y[i]
    end

    if y[end] <= modalpdf
        x[end] = ipdf(y[end])
    else
        x[end] = modalboundary
    end

    x, y
end

## Sampling
# There's a lot of code duplication in this section, but it's difficult to do anything different
# without compromising performance. Be careful when making any changes. Even innocuous looking
# refactors can affect performance. For example, unpacking the Ziggurat and passing it's
# components into another function is faster than passing the Ziggurat object directly. I have
# no idea why.
Base.eltype(::Ziggurat{X}) where {X} = X
Base.eltype(::Type{<:Ziggurat{X}}) where {X} = X

Ytype(::Type{<:Ziggurat{X,Y}}) where {X,Y} = Y
Ytype(::Ziggurat{X,Y}) where {X,Y} = Y

# The where clause is required to force method specialization
@noinline function zigsample_unlikely(parent::P, rng, w, k, y, mb, pdf::F, fb::Nothing, LM, l, x) where {P,F}
    @inbounds begin
        # check density
        yy = (y[l + 1] - y[l]) * rand(rng, eltype(y)) + y[l]
        if yy < pdf(x)
            return x
        end

        # reject sample and retry
        parent(rng, w, k, y, mb, pdf, fb, LM)
    end
end

# The where clause is required to force method specialization
@noinline function zigsample_unlikely(parent::P, rng, w, k, y, mb, pdf::F, fb::FB, LM, l, x) where {P,F,FB}
    @inbounds begin
        if l == 1
            # unbounded tail fallback
            x2 = w[2] * significand_bitmask(eltype(w)) + mb
            return fb(rng, x2)
        end
        # check density
        yy = (y[l + 1] - y[l]) * rand(rng, eltype(y)) + y[l]
        if yy < pdf(x)
            return x
        end

        # reject sample and retry
        parent(rng, w, k, y, mb, pdf, fb, LM)
    end
end

# The where clause is required to force method specialization
@inline function zigsample_floats_masked(rng, w, k, y, mb, pdf::F, fb::FB, LM) where {F,FB}
    r = rand(rng, corresponding_uint(eltype(w)))
    _zigsample_floats_masked(rng, w, k, y, mb, pdf, fb, LM, r)
end

# The where clause is required to force method specialization
@inline function _zigsample_floats_masked(rng, w, k, y, mb, pdf::F, fb::FB, LM, r) where {F,FB}
    l = layer_bits(eltype(w), LM, r) + 1
    u = r >>> shiftbits(eltype(w))
    @inbounds begin
        x = u * w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(zigsample_floats_masked, rng, w, k, y, mb, pdf, fb, LM, l, x)
    end
end

# The where clause is required to force method specialization
@inline function zigsample_floats(rng, w, k, y, mb, pdf::F, fb::FB, LM) where {F,FB}
    r = rand(rng, corresponding_uint(eltype(w)))
    _zigsample_floats(rng, w, k, y, mb, pdf, fb, LM, r)
end

# The where clause is required to force method specialization
@inline function _zigsample_floats(rng, w, k, y, mb, pdf::F, fb::FB, LM, r) where {F,FB}
    l = rand(rng, 1:(length(w) - 1))
    u = r >>> shiftbits(eltype(w))
    @inbounds begin
        x = u * w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(zigsample_floats, rng, w, k, y, mb, pdf, fb, LM, l, x)
    end
end

# The where clause is required to force method specialization
@inline function zigsample_general(rng, w, k, y, mb, pdf::F, fb::FB, LM) where {F,FB}
    l = rand(rng, 1:(length(w) - 1))
    u = rand(rng, eltype(w))
    @inbounds begin
        x = u * w[l] + mb
        if u <= k[l]
            return x
        end
        zigsample_unlikely(zigsample_general, rng, w, k, y, mb, pdf, fb, LM, l, x)
    end
end

function Base.rand(
    rng::AbstractRNG,
    zig_sampler::Random.SamplerTrivial{<:MonotonicZiggurat{X,Y,LM}}
) where {X<:FloatXX,Y,LM}
    z = zig_sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    zigsample_floats_masked(rng, w, k, y, mb, pdf, fb, LM)
end

function Base.rand(
    rng::AbstractRNG,
    zig_sampler::Random.SamplerTrivial{<:MonotonicZiggurat{X,Y,nothing}}
) where {X<:FloatXX,Y}
    z = zig_sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    zigsample_floats(rng, w, k, y, mb, pdf, fb, nothing)
end

function Base.rand(rng::AbstractRNG, zig_sampler::Random.SamplerTrivial{<:MonotonicZiggurat{X,Y}}) where {X,Y}
    z = zig_sampler[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    zigsample_general(rng, w, k, y, mb, pdf, fb, nothing)
end

function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:MonotonicZiggurat{X,Y,LM}}
) where {X<:FloatXX,Y,LM}
    z = s[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    if length(A) < 7 # TODO: Tune this number
        for i in eachindex(A)
            @inbounds A[i] = rand(rng, s)
        end
    else
        T = corresponding_uint(X)
        # UnsafeView is an internal implementation detail of Random.jl
        GC.@preserve A rand!(rng, Random.UnsafeView{T}(pointer(A), length(A)))

        for i in eachindex(A)
            @inbounds r = reinterpret(T, A[i])
            @inbounds A[i] = _zigsample_floats_masked(rng, w, k, y, mb, pdf, fb, LM, r)
        end
    end
    A
end

function Random.rand!(
    rng::Union{TaskLocalRNG,Xoshiro,MersenneTwister},
    A::Array{X},
    s::Random.SamplerTrivial{<:MonotonicZiggurat{X,Y,nothing}}
) where {X<:FloatXX,Y}
    z = s[]
    w = widths(z)
    k = layerratios(z)
    y = heights(z)
    mb = highside(z)
    pdf = density(z)
    fb = fallback(z)
    if length(A) < 7 # TODO: Tune this number
        for i in eachindex(A)
            @inbounds A[i] = rand(rng, s)
        end
    else
        T = corresponding_uint(X)
        # UnsafeView is an internal implementation detail of Random.jl
        GC.@preserve A rand!(rng, Random.UnsafeView{T}(pointer(A), length(A)))

        for i in eachindex(A)
            @inbounds r = reinterpret(T, A[i])
            @inbounds A[i] = _zigsample_floats(rng, w, k, y, mb, pdf, fb, nothing, r)
        end
    end
    A
end

# Fallback for non-FloatXX types is to do nothing (not really a bitmask in that case)
significand_bitmask(::Type{T}) where {T} = oneunit(T)
significand_bitmask(::Type{Float64}) = 0x000fffffffffffff
significand_bitmask(::Type{Float32}) = 0x007fffff
significand_bitmask(::Type{Float16}) = 0x03ff

function fixedbit_fraction(frac)
    if !between(0.0, 1.0, frac)
        error("successive layer width radio is not between zero and one. Got frac=$frac. \
        This means the ziggurat was not constructed correctly.")
    end
    if frac === -zero(frac)
        frac = zero(frac)
    end
    _fixedbit_fraction(frac)
end
# fallback is just a fraction (not really fixedbit in that case)
_fixedbit_fraction(frac) = frac
_fixedbit_fraction(frac::Float64) = reinterpret(Normed{UInt64,52}(frac))
_fixedbit_fraction(frac::Float32) = reinterpret(Normed{UInt32,23}(frac))
_fixedbit_fraction(frac::Float16) = reinterpret(Normed{UInt16,10}(frac))
