using ZigguratTools, Test

include("../utils/Normal.jl")
using .Normal

N = 10
zs = ZigguratSampler(f, finv, F, N, tailsample)

x = sampleziggurat(zs)
@test x isa Float64
@test !isnan(x)
@test !isinf(x)
@test x >= 0

x = symmetricsampleziggurat(zs)
@test x isa Float64
@test !isnan(x)
@test !isinf(x)