plotziggurat(args...) = plotziggurat!(plot(), args...)
plotziggurat!(args...) = plotziggurat!(current(), args...)

function plotziggurat!(p::Plots.Plot, zs::Ziggurat)
    plotziggurat_unboundedsteps!(p, zs)
    plotziggurat_pdf!(p, zs)
end

plotziggurat_pdf!(p::Plots.Plot, zs::Ziggurat) = plot!(p, zs.pdf; color = :blue, lw = 2)

function plotziggurat_unboundedsteps!(p::Plots.Plot, zs::Ziggurat)
    plotziggurat_unboundedsteps!(p, zs.x, zs.y)
end
function plotziggurat_unboundedsteps!(p::Plots.Plot, x::AbstractVector, y::AbstractVector)
    plot!(p, [0, x[1]], [y[1], y[1]]; color = :black, legend = false)
    scatter!(p, [x[1]], [y[1]]; color = :black)
    for i in eachindex(x)[2:end]
        plot!(p, [0, x[i - 1], x[i - 1]], [y[i], y[i], y[i - 1]]; color = :black)
        plot!(p, [x[i], x[i]], [y[i], y[i - 1]]; color = :black, ls = :dash)
        scatter!(p, [x[i]], [y[i]]; color = :black)
    end
    xl = xlims(p)
    yl = ylims(p)
    xlims!(0, 1.5xl[2])
    ylims!(0, yl[2])

    p
end

plotziggurat_boundedsteps(x, y) = plotziggurat_boundedsteps!(plot(), x, y)
plotziggurat_boundedsteps!(x, y) = plotziggurat_boundedsteps!(current(), x, y)
function plotziggurat_boundedsteps!(p::Plots.Plot, x, y)
    plot!(p, [0, x[1]], [y[1], y[1]]; color = :black, legend = false)
    plot!(p, [x[1], x[1]], [0, y[1]]; color = :black)
    for i in eachindex(x)[2:end]
        plot!(p, [0, x[i]], [y[i], y[i]]; color = :black)
        plot!(p, [x[i], x[i]], [y[i - 1], y[i]]; color = :black)
        plot!(p, [x[i], x[i]], [0, y[i - 1]]; color = :black, ls = :dash)
    end

    xl = xlims(p)
    yl = ylims(p)
    xlims!(0, 1.5xl[2])
    ylims!(0, yl[2])

    p
end
