abstract type MonotonicZiggurat{Mask,Shift,X,Y} <: Ziggurat{X} end

struct BoundedZiggurat{Mask,Shift,X,Y,K,F} <: MonotonicZiggurat{Mask,Shift,X,Y}
    w::Vector{X}
    k::Vector{K}
    y::Vector{Y}
    pdf::F
    modalboundary::X
    function BoundedZiggurat{M,S}(w, k, y, f, mb) where {M,S}
        new{M,S,eltype(w),eltype(y),eltype(k),typeof(f)}(w, k, y, f, mb)
    end
end

struct UnboundedZiggurat{Mask,Shift,X,Y,K,F,FB} <: MonotonicZiggurat{Mask,Shift,X,Y}
    w::Vector{X}
    k::Vector{K}
    y::Vector{Y}
    pdf::F
    modalboundary::X
    fallback::FB
    function UnboundedZiggurat{M,S}(w, k, y, f, mb, fb) where {M,S}
        new{M,S,eltype(w),eltype(y),eltype(k),typeof(f),typeof(fb)}(w, k, y, f, mb, fb)
    end
end

corresponding_uint(::Type{Float64}) = UInt64
corresponding_uint(::Type{Float32}) = UInt32
corresponding_uint(::Type{Float16}) = UInt16

const FloatXX = Union{Float64,Float32,Float16}

mask_shift(::Type, N) = nothing, nothing
function mask_shift(X::Type{<:FloatXX}, N)
    maxpower = 8sizeof(X) - Base.significand_bits(X)
    if N ∈ (2^m for m in 0:maxpower)
        power = Int64(log2(N))
        shift = 8sizeof(X) - power
        mask = corresponding_uint(X)(2^power - 1) << shift
        return mask, shift
    else
        return nothing, nothing
    end
end

"""
Provide error messages for unexpected or invalid results.
"""
struct PDFWrap{F,X}
    f::F
    mb::X
    am::X
end

function (pdf::PDFWrap)(x)
    if !between(pdf.mb, pdf.am, x)
        error("Unexpected Error: attempted to evaluate pdf at x = $x outside the domain.")
    end

    result = pdf.f(x)

    if result < 0
        error("pdf($x) is negative. Check the definition of your pdf. It must \
        return positive numbers everywhere on its domain, including the end points of \
        the domain.")
    elseif isnan(result)
        error("pdf($x) is NaN. Check the definition of your pdf. It must \
        return positive numbers everywhere on its domain, including the end points of \
        the domain.")
    elseif isinf(result)
        error("pdf($x) is infinite. The Marsaglia & Tsang Ziggurat \
        algorithm requires a finite pdf.")
    end

    return result
end

"""
Force the ipdf to respect the domain and codomain.
"""
mutable struct IPDFWrap{F,X,Y}
    const ipdf::F
    const mb::X
    const am::X
    const fmb::Y
    const fam::Y
    clamplog::LogLevel
end

IPDFWrap(ipdf, mb, am, fmb, fam, ll = Debug) = IPDFWrap(ipdf, mb, am, fmb, fam, ll)

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
        if result > ub
            @logmsg ipdf.clamplog "inverse pdf has returned a value outside the domain. The result was \
            clamped to be within the domain, but this may indicate an error. Got \
            ipdf($y) = $result, used ipdf($y) = $ub instead."
            return ub
        elseif result < lb
            @logmsg ipdf.clamplog "inverse pdf has returned a value outside the domain. The result was \
            clamped to be within the domain, but this may indicate an error. Got \
            ipdf($y) = $result, used ipdf($y) = $lb instead."
            return lb
        else
            return result
        end
    end
end

wrapperloglevel(unwrapped, ll) = nothing
wrapperloglevel(wipdf::IPDFWrap, ll) = wipdf.clamplog = ll

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

# TODO: This needs documentation.
PDFWrap(f::NoWrap, args...) = f.f
IPDFWrap(f::NoWrap, args...) = f.f

"""
    monotonic_ziggurat(pdf, domain, [N]; [ipdf, tailarea, fallback_generator, ...])

Constructs a high-performance sampler for a univariate probability distribution defined by a
monotonic probability density function (`pdf`). The pdf must not diverge to infinity anywhere
on the domain, including at the endpoints, but may otherwise be arbitrary — including
discontinuous functions. Optional components such as `ipdf`, `tailarea`, and `fallback_generator`
can be customized to avoid slow numerical methods or to support challenging pdfs. Note that
customizing `ipdf` and `tailarea` will only affect the speed and accuracy of the ziggurat
construction, they are not used in the sampling algorithm. Conversely, customizing the
`fallback_generator` will only affect the speed and accuracy of the samples coming from deep
in the tail of the distribution, it is not used in the ziggurat construction algorithm.

# Arguments
 - `pdf`: The probability density function of the desired distribution. It must be monotonic
    and must not diverge to infinity anywhere on the `domain`, including the endpoints. It
    does not need to be normalized, but the other arguments need to be normalized the same
    way. E.g., `ipdf` should be the inverse of the given `pdf` function and therefore must
    have the same normalization.

 - `domain`: The domain of the pdf. `domain` may be any collection of numbers, but only its
    extrema will be used as the boundaries of the domain. The values will be passed into
    `float()` and promoted to a common type. This common type will be the the type produced
    by sampling from the resulting ziggurat. The domain may be unbounded, for example,
    `(0, Inf)`.

 - `N`: (Optional) The number of layers in the ziggurat. If `N` is a power of two and the
    domain is Float64 with N <= 4096, Float32 with N <= 512, or Float16 with N <= 64, then
    an optimized sampling algorithm is used, improving sampling performance. Normally, `N`
    defaults to 256, but for Float16 domains it defaults to 64.

 - `ipdf`: (Optional) This function is the inverse of the given `pdf` argument. It's used in
    the ziggurat construction algorithm, but not during sampling, so it may affect the
    performance of the initial call to `z = monotonic_ziggurat(...)`, but not subsequent
    calls to `rand(z)`. If no ipdf is provided then the inverse will be calculated numerically
    using a root finding algorithm. See also [`inversepdf`](@ref).

 - `tailarea`: (Optional) A function that computes the area of the tail of the `pdf` when the
    domain is unbounded. If the domain has a finite length, this argument is ignored. The
    `tailarea` argument takes precidence over the `cdf` and `ccdf` arguments. If neither
    `tailarea`, nor an appropriate (c)cdf argument is provided then a numerical integral
    will be computed.

 - `cdf`: (Optional) The cumulative distribution function. It should be normalized the same
    way the pdf is normalized. Superceeded by `tailarea`, and ignored when the domain is
    bounded.

 - `ccdf`: (Optional) The complementary cumulative distribution function. It should be
    normalized the same way the pdf is normalized. Superceeded by `tailarea`, and ignored
    when the domain is bounded.

 - `fallback`: (Optional) A function that takes arguments (rng, a) and produces a random
    variate in the tail beyond `a` (below `a` for increasing distributions, or above `a`
    for decreasing distributions). This is only used on unbounded domains, it is ignored on
    bounded domians. If no fallback is provided, then inverse transform sampling will be used
    by numerically inverting the `tailarea` (which may itself may be a numerical estimate).
    Overriding the fallback can greatly improve the worst-case performance of the sampling
    algorithm, but since the fallback is only called rarely, it will not have a large influence
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
"""
function monotonic_ziggurat(
    pdf,
    domain,
    N = nothing;
    ipdf = inversepdf(pdf, domain),
    tailarea = nothing,
    cdf = nothing,
    ccdf = nothing,
    fallback_generator = nothing
)
    domain = regularize(domain)

    if N === nothing
        N = eltype(domain) == Float16 ? 64 : 256
    end

    a, b = extrema(domain)
    if isinf(a) || isinf(b)
        tailarea = _choose_tailarea_func(pdf, domain, tailarea, cdf, ccdf)
        UnboundedZiggurat(pdf, domain, N; ipdf, tailarea, fallback_generator)
    else
        BoundedZiggurat(pdf, domain, N; ipdf)
    end
end

function BoundedZiggurat(pdf, domain, N; ipdf = nothing)
    domain = extrema(regularize(domain))

    _check_arguments(N, domain)
    modalboundary, argminboundary = _identify_mode(pdf, domain)

    wpdf = PDFWrap(pdf, modalboundary, argminboundary)

    if ipdf === nothing
        ipdf = inversepdf(wpdf, domain)
    end

    wipdf = IPDFWrap(
        ipdf,
        modalboundary,
        argminboundary,
        wpdf(modalboundary),
        wpdf(argminboundary)
    )

    if wpdf(modalboundary) == 0
        error("expected the pdf to be non-zero on at least one boundary.")
    end

    if isinf(domain[1]) || isinf(domain[2])
        error("expected a bounded domain, got domain=$domain.")
    end

    # Build ziggurats using wrapped functions
    x, y = search(N, modalboundary, argminboundary, wpdf, wipdf)

    w = (x .- modalboundary) ./ significand_bitmask(eltype(x))
    k = [
        fixedbit_fraction((x[i + 1] - modalboundary)/(x[i] - modalboundary)) for
        i in 1:(length(x) - 1)
    ]

    # final ziggurat uses the unwrapped function so that there is no effect on
    # sampling performance
    BoundedZiggurat{mask_shift(eltype(domain), N)...}(w, k, y, pdf, modalboundary)
end

function UnboundedZiggurat(
    pdf,
    domain,
    N;
    ipdf = nothing,
    tailarea = nothing,
    fallback_generator = nothing
)
    domain = extrema(regularize(domain))

    _check_arguments(N, domain)
    modalboundary, argminboundary = _identify_mode(pdf, domain)

    wpdf = PDFWrap(pdf, modalboundary, argminboundary)

    if ipdf === nothing
        ipdf = inversepdf(wpdf, domain)
    end

    wipdf = IPDFWrap(
        ipdf,
        modalboundary,
        argminboundary,
        pdf(modalboundary),
        pdf(argminboundary)
    )

    if wpdf(modalboundary) == 0
        error("expected the pdf to be non-zero on at least one boundary.")
    end

    if !isinf(domain[1]) && !isinf(domain[2])
        error("expected an unbounded domain, got domain=$domain.")
    end

    if tailarea === nothing
        # TODO: the tailarea function should come from a 'tool' that has this as default
        # and allows customization e.g. quadgk arguments or other integrators.
        modepdf = wpdf(modalboundary)
        domain_type = typeof(modalboundary)
        range_type = typeof(modepdf)
        error_type = typeof(norm(modepdf))
        segbuf = alloc_segbuf(domain_type, range_type, error_type)

        # TODO: Add error tolerance and check the returned error estimate.
        tailarea = let wpdf=wpdf, segbuf=segbuf, argminboundary=argminboundary
            x -> abs(quadgk(wpdf, x, argminboundary; segbuf)[1])
        end
    end

    # Build ziggurats using wrapped functions
    x, y = search(N, modalboundary, argminboundary, wpdf, wipdf, tailarea)

    if fallback_generator === nothing
        # TODO: fallback_generators should come from a 'tool' with this as default but also allows customization.
        # e.g. pass arguments through inverse to find_zero. Think about reducing layers of indirection.
        x2 = x[2]
        ta = tailarea(x2)
        td = if modalboundary > argminboundary
            (nextfloat(typemin(x2)), x2)
        else
            (x2, prevfloat(typemax(x2)))
        end
        inverse_tailprob = let tailarea = tailarea, ta = ta, td = td
            inversepdf(x -> tailarea(x) / ta, td)
        end

        fallback = rng -> inverse_tailprob(rand(rng, typeof(modalboundary)))
    else
        fallback = fallback_generator(x[2])
    end

    w = (x .- modalboundary) ./ significand_bitmask(eltype(x))
    k = [
        fixedbit_fraction((x[i + 1] - modalboundary)/(x[i] - modalboundary)) for
        i in 1:(length(x) - 1)
    ]

    # final ziggurat uses the unwrapped function so that there is no effect on
    # sampling performance
    UnboundedZiggurat{mask_shift(eltype(domain), N)...}(
        w,
        k,
        y,
        pdf,
        modalboundary,
        fallback
    )
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
            throw(ArgumentError("a ccdf is provided for an increasing pdf, pass cdf or tailarea instead."))
        else
            error("Unreachable error: it should be impossible to throw this error.")
        end
    else
        if ccdf !== nothing
            return ccdf
        elseif cdf !== nothing
            throw(ArgumentError("a cdf is provided for an increasing pdf, pass ccdf or tailarea instead."))
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
    # TODO: Roots.Tracks may change between versions, so we should use an alternative (SciML, or in-house?)
    tracker = Roots.Tracks(eltype(p.x), eltype(p.y))
    ystar = find_zero(ziggurat_residual, y_domain, Roots.Bisection(), p; tracks = tracker)

    # y[N] needs to be either exact or a slightly over
    if tracker.convergence_flag !== :exact_zero
        ystar = tracker.abₛ[end][2]
    end
    # TODO handle non-convergence. Bisection is guarenteed to converge, but not all
    # algorithms are.

    # Set the IPDF wrapper to issue warnings when ipdf produces points outside the domain.
    wrapperloglevel(p.buildargs[3], Warn)

    # TODO: Check if the final build fails
    build!(p.x, p.y, ystar, p.buildargs...)
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
function layerarea(y2, x2, modalboundary, tailarea::Function)
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
Random.eltype(::Type{<:MonotonicZiggurat{M,S,X}}) where {M,S,X} = X
Ytype(::Type{<:MonotonicZiggurat{M,S,X,Y}}) where {M,S,X,Y} = Y
Ytype(::MonotonicZiggurat{M,S,X,Y}) where {M,S,X,Y} = Y

function Base.rand(
    rng::AbstractRNG,
    zig_sampler::Random.SamplerTrivial{<:MonotonicZiggurat}
)
    z = zig_sampler[]
    zigsample(rng, z)
end

# Slower fallback for types that don't have specialized sampling algorithms
function zigsample(rng, z::MonotonicZiggurat)
    l = rand(rng, 1:(length(z.w) - 1))
    u = rand(rng, eltype(z))
    x = u * z.w[l] + z.modalboundary
    if u <= z.k[l] # k and u are just fractions for non-FloatXX types, not integers.
        return x
    end
    slowpath(rng, z, l, x)
end

# Specialized/optimized implementation for floats.
function zigsample(rng, z::MonotonicZiggurat{M,S,F}) where {M,S,F<:FloatXX}
    @inbounds begin
        r = rand(rng, corresponding_uint(F))
        l = random_layer(rng, r, z)
        u = r & significand_bitmask(eltype(z))
        x = u * z.w[l] + z.modalboundary
        if u <= z.k[l]
            return x
        end
        slowpath(rng, z, l, x)
    end
end

# Power of two number of layers is optimized
function random_layer(rng, r, z::MonotonicZiggurat{M,S}) where {M,S}
    ((r & M) >> S) + 1
end

# Slower fallback
function random_layer(rng, r, z::MonotonicZiggurat{nothing,nothing})
    rand(rng, 1:(length(z.w) - 1))
end

@noinline function slowpath(rng, z::UnboundedZiggurat, l, x)
    @inbounds begin
        if l == 1
            return z.fallback(rng)
        end

        y = (z.y[l + 1] - z.y[l]) * rand(rng, Ytype(z)) + z.y[l]
        if y < z.pdf(x)
            return x
        end

        zigsample(rng, z)
    end
end

@noinline function slowpath(rng, z::BoundedZiggurat, l, x)
    @inbounds begin
        y = (z.y[l + 1] - z.y[l]) * rand(rng, Ytype(z)) + z.y[l]
        if y < z.pdf(x)
            return x
        end

        zigsample(rng, z)
    end
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
