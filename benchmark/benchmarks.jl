using ZigguratTools, BenchmarkTools

include("Normal.jl")
using .Normal

const SUITE = BenchmarkGroup()

const Nzig = 256
const Nsamp = 1000
SUITE["sampling"] = BenchmarkGroup()
SUITE["sampling"]["oneside"] = BenchmarkGroup()
SUITE["sampling"]["oneside"]["single"] = 
    @benchmarkable sampleziggurat(zs) setup=(zs=ZigguratSampler(f, finv, F, Nzig, tailsample))
SUITE["sampling"]["oneside"]["multiple"] = 
    @benchmarkable sampleziggurat(zs, Nsamp) setup=(zs=ZigguratSampler(f, finv, F, Nzig, tailsample))
