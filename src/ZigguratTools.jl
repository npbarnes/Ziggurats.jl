module ZigguratTools

using Roots
using Random
using Distributions
using Plots

export UnboundedIncreasingZiggurat, UnboundedDecreasingZiggurat, plotziggurat, plotziggurat!

abstract type AbstractZiggurat{X} <: Sampleable{Univariate, Continuous} end
abstract type AbstractMonotonicZiggurat{X} <: AbstractZiggurat{X} end
Base.eltype(::AbstractZiggurat{X}) where X = X

struct UnboundedDecreasingZiggurat{X, Y, F, XY, T} <: AbstractMonotonicZiggurat{X}
    x::Vector{X}
    y::Vector{Y}
    pdf::F
    layerarea::XY
    tailmap::T

    # TODO: add a fallback for tailsample that uses inverse transform sampling
    # or ApproxFun.jl (Olver & Townsend 2013 algorithm)
    # or one or more of the algorithms in Jalalvand & Charsooghi
    function UnboundedDecreasingZiggurat(pdf, ipdf, ccdf, mode, N, tailmap)
        x,y = searchrightziggurat(pdf, ipdf, ccdf, mode, N)
        A = abs(x[1] - mode) * (y[2] - y[1])
        tm = u -> tailmap(x[1], u)
        new{eltype(x), eltype(y), typeof(pdf), typeof(A), typeof(tm)}(x, y, pdf, A, tm)
    end
end

struct UnboundedIncreasingZiggurat{X, Y, F, XY, T} <: AbstractMonotonicZiggurat{X}
    x::Vector{X}
    y::Vector{Y}
    pdf::F
    layerarea::XY
    tailmap::T

    # TODO: add a fallback for tailsample that uses inverse transform sampling
    # or ApproxFun.jl (Olver & Townsend 2013 algorithm)
    # or one or more of the algorithms in Jalalvand & Charsooghi
    function UnboundedIncreasingZiggurat(pdf, ipdf, cdf, mode, N, tailmap)
        x,y = searchleftziggurat(pdf, ipdf, cdf, mode, N)
        A = abs(x[1] - mode) * (y[2] - y[1])
        tm = u -> tailmap(x[1], u)
        new{eltype(x), eltype(y), typeof(pdf), typeof(A), typeof(tm)}(x, y, pdf, A, tm)
    end
end

plotziggurat(zs) = plotziggurat!(plot(), zs)
plotziggurat!(zs) = plotziggurat!(current(), zs)
function plotziggurat!(p::Plots.Plot, zs)
    x, y = zs.x, zs.y

    p = plot!(p, [0, x[1]], [y[1], y[1]], color=:black, legend=false)
    scatter!(p, [x[1]], [y[1]], color=:black)
    for i in eachindex(x)[2:end]
        plot!(p, [0, x[i-1], x[i-1]], [y[i], y[i], y[i-1]], color=:black)
        plot!(p, [x[i], x[i]], [y[i], y[i-1]], color=:black, ls=:dash)
        scatter!(p, [x[i]], [y[i]], color=:black)
    end
    xl = xlims(p)
    yl = ylims(p)
    xlims!(0,1.5xl[2])
    ylims!(0,yl[2])

    plot!(p, zs.pdf, color=:blue, lw=2)

    p
end

Distributions.mode(z::AbstractMonotonicZiggurat) = z.x[end]
Random.rand(rng::AbstractRNG, z::AbstractMonotonicZiggurat) = sampleziggurat(rng, z)

# TODO: add fallback for F(x) = quadgk(f, -Inf, x)
# TODO: add fallback for f(x) = ForwardDiff.derivative(F, x)
# TODO: add interface for supplying F or 1-F
# TODO: add fallback for 1-F(x) = quadgk(f, x, Inf)
# TODO: add fallback for finv using root finding
# TODO: support ziggurats with < 3 layers
function buildziggurat!(x, y, pdf, ipdf, tailarea, mode, x1, mode_pd=pdf(mode))
    x[1] = x1
    y[1] = pdf(x1)
    A = abs(x[1] - mode) * y[1] + tailarea(x[1])
    y[2] = y[1] + A/abs(x[1] - mode)
    for i in eachindex(x)[2:end-1]
        if y[i] >= mode_pd
            # Building the ziggurat failed
            return nothing 
        end
        x[i] = ipdf(y[i])
        y[i+1] = y[i] + A/abs(x[i] - mode)
    end

    x[end] = mode
    
    x,y
end

searchleftziggurat(pdf, ipdf_left, cdf, mode, N) = searchziggurat((-Inf, mode), pdf, ipdf_left, cdf, mode, N)
searchrightziggurat(pdf, ipdf_right, ccdf, mode, N) = searchziggurat((mode, Inf), pdf, ipdf_right, ccdf, mode, N)

function searchziggurat(domain, pdf, ipdf, tailarea, mode, N)
    # TODO: The search method can be made faster since we know x1 is probably closer to 
    # mode than to Inf or -Inf, and we may be able to use some derivative information.
    max_pd = pdf(mode)
    x = Vector{typeof(mode)}(undef, N)
    y = Vector{typeof(max_pd)}(undef, N)

    "Try to build a ziggurat and return a value that represents error"
    function attemptziggurat!(x1)
        zig = buildziggurat!(x, y, pdf, ipdf, tailarea, mode, x1, max_pd)
        if zig === nothing
            return Inf
        end

        x,y = zig
        y[N] - max_pd
    end

    tracker = Roots.Tracks()
    xstar = find_zero(attemptziggurat!, domain, Bisection(); tracks=tracker)

    # y[N] needs to be either exact or a slightly over
    if tracker.convergence_flag !== :exact_zero
        xstar = tracker.abₛ[end][1]
    end
    #TODO handle non-convergence. Bisection is guarenteed to converge, but not all
    #algorithms are.

    buildziggurat!(x, y, pdf, ipdf, tailarea, mode, xstar, max_pd)
end

function sampleziggurat(rng::AbstractRNG, zs::UnboundedDecreasingZiggurat)
    N = length(zs.x)

    while true
        l = rand(rng, 1:N)
        if l == 1 # Baselayer
            δ = zs.layerarea/zs.y[1]
            x = mode(zs) + δ*rand(rng)
            if x < zs.x[1]
                return x
            else
                return zs.tailmap(rand(rng)) # Can I replace rand() with (x-x1)/(m+δ-x1)?
            end
        else # all other layers
            u = rand(rng)
            x = mode(zs)*u + zs.x[l-1]*(1-u)
            if x < zs.x[l]
                return x
            else
                y = (zs.y[l] - zs.y[l-1])*rand(rng) + zs.y[l-1]
                if y < zs.pdf(x)
                    return x
                end
            end
        end
    end
end

function sampleziggurat(rng::AbstractRNG, zs::UnboundedIncreasingZiggurat)
    N = length(zs.x)

    while true
        l = rand(rng, 1:N)
        if l == 1 # Baselayer
            δ = zs.layerarea/zs.y[1]
            x = mode(zs) - δ*rand(rng)
            if zs.x[1] < x
                return x
            else
                return zs.tailmap(rand(rng)) # Can I replace rand() with (x1-x)/(x1-(m-δ))?
            end
        else # all other layers
            u = rand(rng)
            x = mode(zs)*u + zs.x[l-1]*(1-u)
            if zs.x[l] < x
                return x
            else
                y = (zs.y[l] - zs.y[l-1])*rand(rng) + zs.y[l-1]
                if y < zs.pdf(x)
                    return x
                end
            end
        end
    end
end

end # module
