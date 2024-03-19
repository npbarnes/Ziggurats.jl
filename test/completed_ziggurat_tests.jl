using ZigguratTools, Test

include("Normal.jl")
using .Normal

zig = searchziggurat(f, finv, F, 10)

x,y = xyvalues(zig)
A = layerarea(zig)
@test x[1]*y[1] + (1 - F(x[1])) ≈ A
@test all(x[i] * (y[i+1] - y[i]) ≈ A for i in eachindex(x)[1:end-1])

@test all(y ≈ f(x) for (x,y) in zip(x,y))