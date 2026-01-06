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
    q = T == Float16 ? 1e-8 : 1e-6
    nbins = T == Float16 ? 15 : 50
end

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
    q = T == Float16 ? 1e-8 : 1e-6
    nbins = T == Float16 ? 15 : 50
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
                name = "Cauchy (x>0) $T",
                dist = truncated(Cauchy(); lower = 0),
                f = x -> 2 / (π * (1 + x^2)),
                ipdf = y -> √(1/(π * y/2) - 1),
                tailarea = x -> 1 - 2atan(x)/π,
                fallback = (rng, x) -> tan(atan(x) + rand(rng, T)*(T(π)/2 - atan(x))),
                domain = (T(0.0), T(Inf)),
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
                name = "Uniform $T",
                dist = Uniform(),
                f = x -> T(1),
                ipdf = y -> T(1), # Uniform is decreasing
                tailarea = nothing,
                fallback = nothing,
                domain = (T(0), T(1)),
                constructor = BoundedZiggurat,
                T = T
            ),
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

# Laplace is a bell distribution, but It's separated out for extra tests
const LaplaceTestCases = [
    CompositeTestData(;
        name = "Laplace $T",
        dist = Laplace(),
        f = x -> exp(-abs(x))/2,
        ipdfs = [y -> log(2y), y -> -log(2y)],
        cdf = x -> x < 0 ? exp(x)/2 : 1 - exp(-x)/2,
        ccdf = x -> x < 0 ? 1 - exp(x)/2 : exp(-x)/2,
        left_fallback = (rng, x) -> x + log1p(-rand(rng, T)),
        right_fallback = (rng, x) -> x - log1p(-rand(rng, T)),
        domain = (T(-Inf), T(0.0), T(Inf)),
        T = T
    ) for T in TestTypes
]

const OtherBellTestCases = reduce(
    vcat,
    [
        [
            let cdf = x -> (1 + erf(x/T(√2)))/2, ccdf = x -> erfc(x/T(√2))/2
                CompositeTestData(;
                    name = "Normal $T",
                    dist = Normal(),
                    f = x -> 1/√T(2π) * exp(-x^2/2),
                    ipdfs = [(y -> -√(-2log(√T(2π)*y))), (y -> √(-2log(√T(2π)*y)))],
                    cdf = cdf,
                    ccdf = ccdf,
                    left_fallback = (rng, x) -> T(√2 * erfinv(2cdf(x)*(1-rand(rng)) - 1)),
                    right_fallback = (rng, x) -> T(√2 * erfcinv(2ccdf(x)*(1-rand(rng)))),
                    domain = (T(-Inf), T(0.0), T(Inf)),
                    T = T
                )
            end,
            CompositeTestData(;
                name = "Cauchy $T",
                dist = Cauchy(),
                f = x -> 1 / (π * (1 + x^2)),
                ipdfs = [y -> -√(1/(π * y) - 1), y -> √(1/(π * y) - 1)],
                cdf = x -> atan(x)/π + 1/2,
                ccdf = x -> 1/2 - atan(x)/π,
                left_fallback = (rng, x) -> -tan(atan(-x) + rand(rng, T)*(T(π)/2 - atan(-x))),
                right_fallback = (rng, x) -> tan(atan(x) + rand(rng, T)*(T(π)/2 - atan(x))),
                domain = (T(-Inf), T(0.0), T(Inf)),
                T = T
            )
        ] for T in TestTypes
    ]
)

const BellTestCases = [
    LaplaceTestCases;
    OtherBellTestCases
]

const SymmetricTestCases = reduce(
    vcat,
    [
        [let cdf = x -> (1 + erf(x/T(√2)))/2, ccdf = x -> erfc(x/T(√2))/2
            CompositeTestData(;
                name = "Binormal $T",
                dist = MixtureModel(Normal, [(-2, 1), (2, 1)]),
                f = x -> 1/√T(2π) * (exp(-(x-2)^2/2) + exp(-(x+2)^2/2)),
                ipdfs = nothing,
                cdf = x -> cdf(x-2) + cdf(x+2),
                ccdf = x -> ccdf(x-2) + ccdf(x+2),
                left_fallback = (rng, x) -> begin
                    lprob = cdf(x + 2)
                    rprob = cdf(x - 2)
                    p = lprob / (lprob + rprob)
                    if rand(rng, Bernoulli(p))
                        μ = -2
                    else
                        μ = 2
                    end
                    T(√2 * erfinv(2cdf(x - μ)*(1-rand(rng)) - 1)) + μ
                end,
                right_fallback = (rng, x) -> begin
                    lprob = ccdf(x + 2)
                    rprob = ccdf(x - 2)
                    p = rprob / (lprob + rprob)
                    if rand(rng, Bernoulli(p))
                        μ = 2
                    else
                        μ = -2
                    end
                    T(√2 * erfcinv(2ccdf(x - μ)*(1-rand(rng)))) + μ
                end,
                domain = (T(-Inf), T(-1.9986513460302164), T(0.0), T(1.9986513460302164), T(Inf)),
                T = T
            )
        end] for T in TestTypes
    ]
)

const AsymmetricTestCases = reduce(
    vcat,
    [
        [let k = 3, A = 1/(2^(T(k)/2) * gamma(T(k)/2))
            CompositeTestData(;
                name = "Chisq($k) $T",
                dist = Chisq(k),
                f = x -> isinf(x) ? zero(x) : A * √x * exp(-x/2),
                ipdfs = [y -> -lambertw(-(y/A)^2, 0), y -> -lambertw(-(y/A)^2, -1)],
                cdf = x -> gamma_inc(1.5, x/2)[1],
                ccdf = x -> gamma_inc(1.5, x/2)[2],
                left_fallback = nothing,
                right_fallback = (rng, x0) -> begin
                    if x0 < 2
                        while true
                            x = randn(rng, T)^2 + randn(rng, T)^2 + randn(rng, T)^2
                            if x > x0
                                return x
                            end
                        end
                    end
                    β = 2/(x0 - 1)
                    while true
                        z = β * randexp(rng, T)
                        u = rand(rng, T)
                        if u^2 <= (1 + z) * exp(-z)
                            return x0 * (1 + z)
                        end
                    end
                end,
                domain = (T(0), T(1), T(Inf)),
                T = T
            )
        end] for T in TestTypes
    ]
)

const CompositeTestCases = [
    BellTestCases;
    SymmetricTestCases;
    AsymmetricTestCases
]
