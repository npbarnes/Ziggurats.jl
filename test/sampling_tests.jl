@testset "Sampling Tests" begin
    @testset "Domain type: $T" for T in (Float16, Float32, Float64)
        @testset "Unbounded Ziggurats" begin
            @testset "Normal (x>=0)" begin
                dist = truncated(Normal(); lower = 0)

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> √(-2log(√T(2π)*y))
                tailarea = x -> erfc(x/T(√2))/2
                fallback = (rng, x) -> T(√2 * erfcinv(2tailarea(x)*(1-rand(rng))))
                domain = (T(0), T(Inf))

                # Covers special cases:
                # - N=1: power of two, single layer
                # - N=2: power of two, no middle layers
                # - N=3: not a power of two, with middle layer
                # - N=4: power of two, with middle layers
                # - default rng, and both built in rng's
                # - array generation and scalar generation
                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = UnboundedZiggurat(f, domain, N; ipdf, tailarea, fallback)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, dist; rng, array_generation)
                end
            end

            @testset "Normal (x<=0)" begin
                dist = truncated(Normal(); upper = 0)

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> -√(-2log(√T(2π)*y))
                tailarea = x -> (1 + erf(x/T(√2)))/2
                fallback = (rng, x) -> T(√2 * erfinv(2tailarea(x)*(1-rand(rng,)) - 1))
                domain = (T(-Inf), T(0))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = UnboundedZiggurat(f, domain, N; ipdf, tailarea, fallback)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, dist; rng, array_generation)
                end
            end

            @testset "Exponential" begin
                dist = Exponential()

                f = x -> exp(-x)
                ipdf = y -> -log(y)
                tailarea = x -> exp(-x)
                fallback = (rng, x) -> x - log1p(-rand(rng, T))
                domain = (T(0), T(Inf))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = UnboundedZiggurat(f, domain, N; ipdf, tailarea, fallback)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, dist; rng, array_generation)
                end
            end

            @testset "SteppedExponential" begin
                dist = SteppedExponential(oneunit(T))

                f = Base.Fix1(pdf, dist)
                ta = Base.Fix1(ccdf, dist)
                domain = (T(0), T(Inf))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = UnboundedZiggurat(f, domain, N; tailarea = ta)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, dist; rng, array_generation)
                end
            end

            @testset "TDist (x>0)" begin
                dist = truncated(TDist(1); lower = 0)

                f = x -> (1 + x^2)^-1
                ipdf = y -> √(y^-1 - 1)
                tailarea = x -> (T(π) - 2atan(x))/2
                fallback = (rng, x) -> tan(-((T(π)-2atan(x))*rand(rng, T) - T(π))/2)
                domain = (T(0), T(Inf))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = UnboundedZiggurat(f, domain, N; ipdf, tailarea, fallback)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    q = T == Float16 ? 5e-7 : 1e-6 # Lower confidence in Float16s
                    test_samples(z, dist; q, rng, array_generation)
                end
            end
        end

        @testset "Bounded Ziggurats" begin
            @testset "Truncated Normal (0.5 <= x <= 1)" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> √(-2log(√T(2π)*y))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = BoundedZiggurat(f, (T(0.5), T(1)), N; ipdf)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, truncated(dist; lower = 0.5, upper = 1); rng, array_generation)
                end
            end

            @testset "Truncated Normal (-1 <= x <= -0.5)" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> -√(-2log(√T(2π)*y))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = BoundedZiggurat(f, (T(-1), T(-0.5)), N; ipdf)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, truncated(dist; lower = -1, upper = -0.5); rng, array_generation)
                end
            end

            @testset "Truncated Normal (0.5 <= x <= 10)" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> √(-2log(√T(2π)*y))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = BoundedZiggurat(f, (T(0.5), T(10)), N; ipdf)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, truncated(dist; lower = 0.5, upper = 10); rng, array_generation)
                end
            end

            @testset "Truncated Normal (-10 <= x <= -0.5)" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> -√(-2log(√T(2π)*y))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = BoundedZiggurat(f, (T(-10), T(-0.5)), N; ipdf)
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, truncated(dist; lower = -10, upper = -0.5); rng, array_generation)
                end
            end
        end

        @testset "Composite Ziggurats" begin
            @testset "Normal" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                left_ipdf = y -> -√(-2log(√T(2π)*y))
                right_ipdf = y -> √(-2log(√T(2π)*y))
                cdf = x -> (1 + erf(x/T(√2)))/2
                ccdf = x -> erfc(x/T(√2))/2
                left_fallback = (rng, x) -> T(√2 * erfinv(2cdf(x)*(1-rand(rng,)) - 1))
                right_fallback = (rng, x) -> T(√2 * erfcinv(2ccdf(x)*(1-rand(rng))))
                domain = (T(-Inf), T(0.0), T(Inf))

                @testset for N in [1, 2, 3, 4], rng in [missing, Xoshiro(1234), MersenneTwister(1234)], array_generation in [true, false]
                    z = CompositeZiggurat(
                        f,
                        domain,
                        N;
                        ipdfs = [left_ipdf, right_ipdf],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == T
                    @test eltype(rand(z, 3)) == T
                    test_samples(z, dist; rng, array_generation)
                end
            end
        end
    end
end

nothing
