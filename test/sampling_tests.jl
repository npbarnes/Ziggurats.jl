@kwdef struct MonotonicTestData
    name::Any
    dist::Any
    f::Any
    ipdf::Any
    tailarea::Any
    fallback::Any
    domain::Any
    constructor::Any
    T::Any
    # Float16 doen't have enough fidelity for 50 bins, one bin ends up being ~40 Float16's wide.
    # Rounding and the underlying unform RNG seem to cause some bins to be over represented and
    # some bins to be under represented. Using fewer bins and lower confidence for Float16's to
    # mitigate the problem.
    q = T == Float16 ? 1e-7 : 1e-6
    nbins = T == Float16 ? 25 : 50
end

const TestTypes = (Float64, Float32, Float16)

const UnboundedTestCases = reduce(
    vcat,
    [
        [
            MonotonicTestData(;
                name = "Normal (x>=0) $T",
                dist = truncated(Normal(); lower = 0),
                f = x -> 1/√T(2π) * exp(-x^2/2),
                ipdf = y -> √(-2log(√T(2π)*y)),
                tailarea = x -> erfc(x/T(√2))/2,
                fallback = (rng, x) -> T(√2 * erfcinv(erfc(x/T(√2))*(1-rand(rng)))),
                domain = (T(0), T(Inf)),
                constructor = UnboundedZiggurat,
                T = T
            ),
            MonotonicTestData(;
                name = "Normal (x<=0) $T",
                dist = truncated(Normal(); upper = 0),
                f = x -> 1/√T(2π) * exp(-x^2/2),
                ipdf = y -> -√(-2log(√T(2π)*y)),
                tailarea = x -> (1 + erf(x/T(√2)))/2,
                fallback = (rng, x) -> T(√2 * erfinv((1 + erf(x/T(√2)))*(1-rand(rng)) - 1)),
                domain = (T(-Inf), T(0)),
                constructor = UnboundedZiggurat,
                T = T
            ),
            MonotonicTestData(;
                name = "Exponential $T",
                dist = Exponential(),
                f = x -> exp(-x),
                ipdf = y -> -log(y),
                tailarea = x -> exp(-x),
                fallback = (rng, x) -> x - log1p(-rand(rng, T)),
                domain = (T(0), T(Inf)),
                constructor = UnboundedZiggurat,
                T = T
            ),
            let dist = SteppedExponential(oneunit(T))
                MonotonicTestData(;
                    name = "SteppedExponential $T",
                    dist = dist,
                    f = Base.Fix1(pdf, dist),
                    ipdf = ipdf_SteppedExponential(dist),
                    tailarea = Base.Fix1(ccdf, dist),
                    fallback = SteppedExpFallback(dist),
                    domain = (T(0), T(Inf)),
                    constructor = UnboundedZiggurat,
                    T = T
                )
            end,
            MonotonicTestData(;
                name = "TDist (x>0) $T",
                dist = truncated(TDist(1); lower = 0),
                f = x -> (1 + x^2)^-1,
                ipdf = y -> √(y^-1 - 1),
                # tailarea and fallback are computed in Float64 and converted for more accuracy.
                # Otherwise sampling tests fail in Float16.
                tailarea = x -> T((π - 2atan(x))/2),
                fallback = (rng, x) -> T(tan(-((π-2atan(x))*rand(rng) - π)/2)),
                domain = (T(0), T(Inf)),
                constructor = UnboundedZiggurat,
                T = T
            )
        ] for T in TestTypes
    ]
)

const BoundedTestCases = reduce(
    vcat,
    [
        [
            MonotonicTestData(;
                name = "Truncated Normal (0.5 <= x <= 1) $T",
                dist = truncated(Normal(); lower = 0.5, upper = 1),
                f = x -> 1/√T(2π) * exp(-x^2/2),
                ipdf = y -> √(-2log(√T(2π)*y)),
                tailarea = nothing,
                fallback = nothing,
                domain = (T(0.5), T(1)),
                constructor = BoundedZiggurat,
                T = T
            ),
            MonotonicTestData(;
                name = "Truncated Normal (-1 <= x <= -0.5) $T",
                dist = truncated(Normal(); lower = -1, upper = -0.5),
                f = x -> 1/√T(2π) * exp(-x^2/2),
                ipdf = y -> -√(-2log(√T(2π)*y)),
                tailarea = nothing,
                fallback = nothing,
                domain = (T(-1), T(-0.5)),
                constructor = BoundedZiggurat,
                T = T
            ),
            MonotonicTestData(;
                name = "Truncated Normal (0.5 <= x <= 10) $T",
                dist = truncated(Normal(); lower = 0.5, upper = 10),
                f = x -> 1/√T(2π) * exp(-x^2/2),
                ipdf = y -> √(-2log(√T(2π)*y)),
                tailarea = nothing,
                fallback = nothing,
                domain = (T(0.5), T(10)),
                constructor = BoundedZiggurat,
                T = T
            ),
            MonotonicTestData(;
                name = "Truncated Normal (-10 <= x <= -0.5) $T",
                dist = truncated(Normal(); lower = -10, upper = -0.5),
                f = x -> 1/√T(2π) * exp(-x^2/2),
                ipdf = y -> -√(-2log(√T(2π)*y)),
                tailarea = nothing,
                fallback = nothing,
                domain = (T(-10), T(-0.5)),
                constructor = BoundedZiggurat,
                T = T
            )
        ] for T in TestTypes
    ]
)

const MonotonicTestCases = [
    UnboundedTestCases;
    BoundedTestCases
]

@kwdef struct CompositeTestData
    name::Any
    dist::Any
    f::Any
    ipdfs::Any
    cdf::Any
    ccdf::Any
    left_fallback::Any
    right_fallback::Any
    domain::Any
    T::Any
    q = 1e-6
end

const SymmetricTestCases = reduce(
    vcat,
    [
        [let cdf = x -> (1 + erf(x/T(√2)))/2, ccdf = x -> erfc(x/T(√2))/2
            CompositeTestData(;
                name = "Normal $T",
                dist = Normal(),
                f = x -> 1/√T(2π) * exp(-x^2/2),
                ipdfs = [(y -> -√(-2log(√T(2π)*y))), (y -> √(-2log(√T(2π)*y)))],
                cdf = cdf,
                ccdf = ccdf,
                left_fallback = (rng, x) -> T(√2 * erfinv(2cdf(x)*(1-rand(rng,)) - 1)),
                right_fallback = (rng, x) -> T(√2 * erfcinv(2ccdf(x)*(1-rand(rng)))),
                domain = (T(-Inf), T(0.0), T(Inf)),
                T = T
            )
        end] for T in TestTypes
    ]
)

const CompositeTestCases = [SymmetricTestCases;]

function monotonic_sampling_tests(td::MonotonicTestData)
    @testset for N in [1, 2, 3, 4]
        @testset "rng isa $(typeof(_rng))" for _rng in [missing, Xoshiro(1234), MersenneTwister(1234)]
            @testset for array_generation in [true, false]
                test_fallbacks = td.fallback === nothing ? [nothing] : [nothing, td.fallback]
                @testset "fallback provided: $(fallback !== nothing)" for fallback in test_fallbacks
                    z = monotonic_ziggurat(td.f, td.domain, N; td.ipdf, td.tailarea, fallback)

                    @test z isa td.constructor
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == td.T
                    @test eltype(rand(z, 3)) == td.T

                    test_samples(z, td.dist; td.nbins, td.q, rng = deepcopy(_rng), array_generation)
                end
            end
        end
    end
end

function composite_sampling_tests(td::CompositeTestData)
    @testset for N in [1, 2, 3, 4]
        @testset "rng isa $(typeof(_rng))" for _rng in [missing, Xoshiro(1234), MersenneTwister(1234)]
            @testset for array_generation in [true, false]
                test_leftfallbacks = td.left_fallback === nothing ? [nothing] : [nothing, td.left_fallback]
                test_rightfallbacks = td.right_fallback === nothing ? [nothing] : [nothing, td.right_fallback]
                @testset for left_fallback in test_leftfallbacks, right_fallback in test_rightfallbacks
                    z = CompositeZiggurat(td.f, td.domain, N; td.ipdfs, td.cdf, td.ccdf, left_fallback, right_fallback)

                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == td.T
                    @test eltype(rand(z, 3)) == td.T

                    test_samples(z, td.dist; td.q, rng = deepcopy(_rng), array_generation)
                end
            end
        end
    end
end

@testset "Sampling Tests" begin
    @testset "$(td.name)" for td in MonotonicTestCases
        monotonic_sampling_tests(td)
    end
    @testset "$(td.name)" for td in CompositeTestCases
        composite_sampling_tests(td)
    end
end

nothing
