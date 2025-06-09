abstract type MonotonicZiggurat{X} end

struct BoundedZiggurat{X,Y,F} <: MonotonicZiggurat{X}
    x::Vector{X}
    y::Vector{Y}
    pdf::F
    modalboundary::X
end

struct UnboundedZiggurat{X,Y,F,FB} <: MonotonicZiggurat{X}
    x::Vector{X}
    y::Vector{Y}
    pdf::F
    modalboundary::X
    fallback::FB
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

function monotonic_ziggurat(
    pdf,
    domain,
    N = 256;
    ipdf = inverse(pdf, domain),
    tailarea = nothing,
    cdf = nothing,
    ccdf = nothing,
    fallback_generator = nothing
)
    if isinf(domain[1]) || isinf(domain[2])
        tailarea = _choose_tailarea_func(pdf, domain, tailarea, cdf, ccdf)
        UnboundedZiggurat(pdf, domain, N; ipdf, tailarea, fallback_generator)
    else
        BoundedZiggurat(pdf, domain, N; ipdf)
    end
end

function BoundedZiggurat(pdf, domain, N; ipdf = inverse(pdf, domain))
    BoundedZiggurat(pdf, domain, N, ipdf)
end

function BoundedZiggurat(pdf, domain, N, ipdf)
    domain = extrema(regularize_domain(domain))

    _check_arguments(N, domain)
    modalboundary, argminboundary = _identify_mode(pdf, domain)

    if pdf(modalboundary) == 0
        error("expected the pdf to be non-zero on at least one boundary.")
    end

    if isinf(domain[1]) || isinf(domain[2])
        error("expected a bounded domain, got domain=$domain.")
    end

    x, y = search(N, modalboundary, argminboundary, pdf, ipdf)

    BoundedZiggurat(x, y, pdf, modalboundary)
end

function UnboundedZiggurat(
    pdf,
    domain,
    N;
    ipdf = inverse(pdf, domain),
    tailarea = nothing,
    fallback_generator = nothing
)
    UnboundedZiggurat(pdf, N, domain, ipdf, tailarea, fallback_generator)
end

function UnboundedZiggurat(pdf, N, domain, ipdf, tailarea, fallback_generator)
    domain = extrema(regularize_domain(domain))

    _check_arguments(N, domain)
    modalboundary, argminboundary = _identify_mode(pdf, domain)

    if pdf(modalboundary) == 0
        error("expected the pdf to be non-zero on at least one boundary.")
    end

    if !isinf(domain[1]) && !isinf(domain[2])
        error("expected an unbounded domain, got domain=$domain.")
    end

    if tailarea === nothing
        # TODO: the tailarea function should come from a 'tool' that has this as default
        # and allows customization e.g. quadgk arguments or other integrators.
        modepdf = pdf(modalboundary)
        domain_type = typeof(modalboundary)
        range_type = typeof(modepdf)
        error_type = typeof(norm(modepdf))
        segbuf = alloc_segbuf(domain_type, range_type, error_type)

        # TODO: Add error tolerance and check the returned error estimate.
        tailarea = let pdf=pdf, segbuf=segbuf, argminboundary=argminboundary
            x -> abs(quadgk(pdf, x, argminboundary; segbuf)[1])
        end
    end

    x, y = search(N, modalboundary, argminboundary, pdf, ipdf, tailarea)

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
            inverse(x -> tailarea(x) / ta, td)
        end

        fallback = rng -> inverse_tailprob(rand(rng))
    else
        fallback = fallback_generator(x[2])
    end

    UnboundedZiggurat(x, y, pdf, modalboundary, fallback)
end

function _check_arguments(N, domain)
    if N < 1
        throw(DomainError(N, "N must be a positive integer, got N=$N."))
    end

    # Check if the domain is well formed and appropriate for a monotonic
    # distribution. I.e. d[1] < d[2], and at most one of d[1] and d[2] are
    # infinite.
    if domain[1] == domain[2]
        error("empty domains are not allowed, got domain=$domain.")
    end
    if domain[1] > domain[2]
        error("malformed domain. domain[1] must be less than domain[2], got domain=$domain.")
    end
    if isinf(domain[1]) && isinf(domain[2])
        error("a domain of (-Inf, Inf) is impossible for a monotonic distribution.")
    end

    return nothing
end

function _identify_mode(pdf, domain)
    # Return the modalboundary (mb) and argminboundary (am).
    # Assume that the domain is well formed and appropriate for a monotonic
    # distribution. I.e. d[1] < d[2], and at most one of d[1] and d[2] are
    # infinite.
    if isinf(domain[1])
        mb = domain[2]
        am = domain[1]
    elseif isinf(domain[2])
        mb = domain[1]
        am = domain[2]
    else
        boundarypdf = pdf.(domain)
        if boundarypdf[1] == boundarypdf[2]
            error("pdf must be monotonic and non-constant.")
        else
            mb = domain[argmax(boundarypdf)]
            am = domain[argmin(boundarypdf)]
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
    tracker = Roots.Tracks()
    ystar = find_zero(ziggurat_residual, y_domain, Bisection(), p; tracks = tracker)

    # y[N] needs to be either exact or a slightly over
    if tracker.convergence_flag !== :exact_zero
        ystar = tracker.abₛ[end][2]
    end
    # TODO handle non-convergence. Bisection is guarenteed to converge, but not all
    # algorithms are.

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
Random.eltype(::Type{<:MonotonicZiggurat{X}}) where {X} = X

function Base.rand(
    rng::AbstractRNG,
    zig_sampler::Random.SamplerTrivial{<:MonotonicZiggurat}
)
    z = zig_sampler[]
    N = length(z.x) - 1 # number of layers

    while true
        l = rand(rng, 1:N)
        x = (z.modalboundary - z.x[l]) * rand(rng) + z.x[l]

        if between(z.x[l + 1], z.modalboundary, x)
            return x
        else
            sp = slowpath(rng, z, l, x)
            if sp !== nothing
                return sp
            end
        end
    end
end

function simple_rejection(rng, z, l, x)
    y = (z.y[l + 1] - z.y[l]) * rand(rng) + z.y[l]
    if y < z.pdf(x)
        return x
    end

    nothing
end

function slowpath(rng, z::UnboundedZiggurat, l, x)
    if l == 1
        return z.fallback(rng)
    end

    simple_rejection(rng, z, l, x)
end

slowpath(rng, z::BoundedZiggurat, l, x) = simple_rejection(rng, z, l, x)
