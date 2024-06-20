using ZigguratTools, BenchmarkTools

using Normal

suite = BenchmarkGroup()

zigsizes = [128, 256]
sampsizes = [1, 10, 100, 1000]

for Nzig in zigsizes
    zs = ZigguratSampler(f, finv, F, Nzig, tailsample)

    suite["ziggurat"][("Nzig", Nzig)] = 
        @benchmarkable ZigguratSampler($f, $finv, $F, $Nzig, $tailsample)

    suite["sampling"]["oneside"]["single"][("Nzig", Nzig)] =
        @benchmarkable sampleziggurat($zs)
    suite["sampling"]["symmetric"]["single"][("Nzig", Nzig)] =
        @benchmarkable symmetricsampleziggurat($zs)
    
    for Nsamp in sampsizes
        suite["sampling"]["oneside"]["multiple"][("Nzig", Nzig)][("Nsamp", Nsamp)] = 
            @benchmarkable sampleziggurat($zs, $Nsamp)
        suite["sampling"]["symmetric"]["multiple"][("Nzig", Nzig)][("Nsamp", Nsamp)] = 
            @benchmarkable symmetricsampleziggurat($zs, $Nsamp)
    end
end

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "tune.json")

if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals)
else
    tune!(suite)
    BenchmarkTools.save(paramspath, params(suite))
end

results = run(suite; verbose=true)
resultspath = joinpath(dirname(@__FILE__), "result.json")
BenchmarkTools.save(resultspath, results)
results
