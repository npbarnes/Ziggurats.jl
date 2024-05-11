using ZigguratTools, Test

include("../utils/Normal.jl")
using .Normal

Nzig = 10
Nsamp = 5

zs = ZigguratSampler(f, finv, F, Nzig, tailsample)

x = sampleziggurat(zs)
@test x isa Float64
@test !isnan(x)
@test !isinf(x)
@test x >= 0

xs = sampleziggurat(zs, Nsamp)
@test length(xs) == Nsamp
@test xs isa Vector{Float64}
@test all(!isnan, xs)
@test all(!isinf, xs)
@test all(>=(0), xs)

out = Vector{Float64}(undef, Nsamp)
xs = sampleziggurat!(out, zs)
@test xs === out
@test length(xs) == Nsamp
@test xs isa Vector{Float64}
@test all(!isnan, xs)
@test all(!isinf, xs)
@test all(>=(0), xs)

x = symmetricsampleziggurat(zs)
@test x isa Float64
@test !isnan(x)
@test !isinf(x)

xs = symmetricsampleziggurat(zs, Nsamp)
@test length(xs) == Nsamp
@test xs isa Vector{Float64}
@test all(!isnan, xs)
@test all(!isinf, xs)

out = Vector{Float64}(undef, Nsamp)
xs = symmetricsampleziggurat!(out, zs)
@test xs === out
@test length(xs) == Nsamp
@test xs isa Vector{Float64}
@test all(!isnan, xs)
@test all(!isinf, xs)
