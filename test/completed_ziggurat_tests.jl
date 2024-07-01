using ZigguratTools, Test

include("./Normal.jl")

N = 10
zs = ZigguratSampler(f, finv, F, N, tailsample)

x,y = xyvalues(zs)
A = layerarea(zs)
f0 = f(0)

# Number of layers
@test length(x) == length(y) == N 

# y = f(x)
@test all(y ≈ f(x) for (x,y) in zip(x,y)) 

# Base layer area
@test x[1]*y[1] + (1 - F(x[1])) ≈ A 

# Other layer areas
@test all(x[i] * (y[i+1] - y[i]) ≈ A for i in eachindex(x)[1:end-1]) 

# top layer is close to and greater than f(0)
@test y[end] ≈ f0 && y[end] >= f0