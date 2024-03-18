module ZigguratTools

using Roots
using Plots

export buildziggurat, buildziggurat!
export searchziggurat
export sampleziggurat
export plotziggurat
export Ziggurat, xvalues, yvalues, xyvalues, layerarea

struct Ziggurat{X,Y,XY}
    x::Vector{X}
    y::Vector{Y}
    A::XY
end

xvalues(z::Ziggurat) = z.x
yvalues(z::Ziggurat) = z.y
xyvalues(z::Ziggurat) = (z.x, z.y)
layerarea(z::Ziggurat) = z.A

function plotziggurat(zig::Ziggurat)
    x,y = xyvalues(zig)

    p = plot([0, x[1]], [y[1], y[1]], color=:black, legend=false)
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
    p
end
function plotziggurat(zig::Ziggurat, f)
    p = plotziggurat(zig)
    xl = xlims(p)
    x = range(xl[1],xl[2],length=1000)
    plot!(p, x, f.(x), color=:blue, lw=2)
end

function buildziggurat(f, finv, F, N, x1; f0=f(zero(x1)))
    x = Vector{typeof(x1)}(undef, N)
    y = Vector{typeof(f0)}(undef, N)
    buildziggurat!(x, y, f, finv, F, x1; f0)
end

"""
    buildziggurat!(x, y, f, finv, F, N::Integer, x1; f0=f(zero(eltype(x))))

Build a table of ziggurat x and y values. Raise and error if the ziggurat is invalid.
"""
function buildziggurat!(x, y, f, finv, F, x1; f0=f(zero(x1)))
    # Check arguments
    if !isreal(x1) || x1 <= 0
        throw(ArgumentError("x1 must be a real positive number."))
    end
    if length(x) != length(y)
        throw(ArgumentError("the lengths of x and y must be the same"))
    end

    # Do the actual work
    zig = _buildziggurat!(x, y, f, finv, F, x1; f0)
    
    # Sanity checks
    if zig === nothing
        error("failed to build a covering ziggurat. The starting value, x1=$x1, is too small for N=$N.")
    end

    x,y = xyvalues(zig)
    if y[end] < f0
        error("failed to buld a covering ziggurat. The starting value, x1=$x1, is too large for N=$N.")
    elseif y[end] ≉ f0 # i.e y[end] >> f0
        error("failed to build an accurate ziggurat. y[N] = $(y[N]) is not close to f(0) = $f0.\
        The starting value, x1=$x1, is too small for N=$N.")
    end

    zig
end

function _buildziggurat!(x, y, f, finv, F, x1; f0)
    N = length(x) # length(x) == length(y)
    x[1] = x1
    y[1] = f(x1)
    A = x[1] * y[1] + (1 - F(x[1]))

    y[2] = y[1] + A/x[1]
    for i in 2:N-1
        if y[i] >= f0
            # Building the ziggurat failed
            return nothing 
        end
        x[i] = finv(y[i])
        y[i+1] = y[i] + A/x[i]
    end

    # The sampling algorithm requires x[N] = 0.
    x[N] = zero(eltype(x))
    
    Ziggurat(x,y,A)
end


"""
Build a ziggurat for f when x1 is unknown by searching for the correct value.
Returns the completed ziggurat table of x's and y's.
"""
function searchziggurat(f, finv, F, N, method=Bisection(); xtype=Float64)
    f0 = f(zero(xtype))
    x = Vector{xtype}(undef, N)
    y = Vector{typeof(f0)}(undef, N)

    "Try to build a ziggurat and return a value that represents error"
    function attemptziggurat!(x1)
        zig = _buildziggurat!(x, y, f, finv, F, x1; f0)
        if zig === nothing
            return Inf
        end

        yvalues(zig)[N] - f0
    end

    tracker = Roots.Tracks()
    xstar = find_zero(attemptziggurat!, (0, Inf), method; tracks=tracker)

    # y[N] needs to be either exact or a slight over estimate
    if tracker.convergence_flag !== :exact_zero
        xstar = tracker.abₛ[end][1]
    end
    #TODO handle non-convergence.

    buildziggurat!(x, y, f, finv, F, xstar; f0)
end

function sampleziggurat!(out, zig)
    for i in eachindex(out)
        out[i] = sampleziggurat(zig)
    end
    out
end

function sampleziggurat(zig, N)
    out = Vector{Float64}(undef, N)
    sampleziggurat!(out, zig)
end

function sampleziggurat(zig)
    N = length(zig.x)

    while true
        l = rand(1:N)
        if l == 1 # Baselayer
            x0 = zig.A/zig.y[1]
            if x0 < zig.x[1]
                return x0
            else
                return tailsample(zig.x[1])
            end
        else # all other layers
            x = zig.x[l-1]*rand()
            if x < zig.x[l]
                return x
            else
                y = (zig.y[l] - zig.y[l-1])*rand() + zig.y[l-1]
                if y < f(x)
                    return x
                end
            end
        end
    end
end

end # module ZigguratTools
