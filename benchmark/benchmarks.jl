using ZigguratTools, BenchmarkTools

include("../utils/Normal.jl")
using .Normal

const SUITE = BenchmarkGroup()

zigsizes = [128, 256]
sampsizes = [1, 10, 100, 1000]

SUITE["sampling"] = BenchmarkGroup()
SUITE["ziggurat"] = BenchmarkGroup()
SUITE["sampling"]["oneside"] = BenchmarkGroup()
SUITE["sampling"]["symmetric"] = BenchmarkGroup()
SUITE["sampling"]["oneside"]["single"] = BenchmarkGroup()
SUITE["sampling"]["symmetric"]["single"] = BenchmarkGroup()

for Nzig in zigsizes
    zs = ZigguratSampler(f, finv, F, Nzig, tailsample)

    SUITE["ziggurat"][("Nzig", Nzig)] = 
        @benchmarkable ZigguratSampler($f, $finv, $F, $Nzig, $tailsample)

    SUITE["sampling"]["oneside"]["single"][("Nzig", Nzig)] =
        @benchmarkable sampleziggurat($zs)
    SUITE["sampling"]["symmetric"]["single"][("Nzig", Nzig)] =
        @benchmarkable symmetricsampleziggurat($zs)

    SUITE["sampling"]["oneside"]["multiple"][("Nzig", Nzig)] = BenchmarkGroup()
    SUITE["sampling"]["symmetric"]["multiple"][("Nzig", Nzig)] = BenchmarkGroup()
    
    for Nsamp in sampsizes
        SUITE["sampling"]["oneside"]["multiple"][("Nzig", Nzig)][("Nsamp", Nsamp)] = 
            @benchmarkable sampleziggurat($zs, $Nsamp)
        SUITE["sampling"]["symmetric"]["multiple"][("Nzig", Nzig)][("Nsamp", Nsamp)] = 
            @benchmarkable symmetricsampleziggurat($zs, $Nsamp)
    end
end
