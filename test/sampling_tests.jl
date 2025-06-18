@testset "Sampling Tests" begin
    @testset "Domain type: $T" for T in (Float16, Float32, Float64)
        @testset "Unbounded Ziggurats" begin
            @testset "Normal (x>=0)" begin
                dist = truncated(Normal(); lower = 0)

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> √(-2log(√T(2π)*y))
                tailarea = x -> erfc(x/T(√2))/2
                fallback_generator = x -> rng -> T(√2 * erfcinv(2tailarea(x)*(1-rand(rng))))
                domain = (T(0), T(Inf))

                # Covers special cases:
                # - N=1: power of two, single layer
                # - N=2: power of two, no middle layers
                # - N=3: not a power of two, with middle layer
                # - N=4: power of two, with middle layers
                @testset for N in [1, 2, 3, 4]
                    z = UnboundedZiggurat(f, domain, N; ipdf, tailarea, fallback_generator)
                    @test typeof(rand(z)) == eltype(z) == ZigguratTools.Ytype(z) == T
                    test_samples(z, dist)
                end
            end

            @testset "Normal (x<=0)" begin
                dist = truncated(Normal(); upper = 0)

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> -√(-2log(√T(2π)*y))
                tailarea = x -> (1 + erf(x/T(√2)))/2
                fallback_generator =
                    x -> rng -> T(√2 * erfinv(2tailarea(x)*(1-rand(rng,)) - 1))
                domain = (T(-Inf), T(0))

                @testset for N in [1, 2, 3, 4]
                    z = UnboundedZiggurat(f, domain, N; ipdf, tailarea, fallback_generator)
                    @test typeof(rand(z)) == eltype(z) == ZigguratTools.Ytype(z) == T
                    test_samples(z, dist)
                end
            end

            @testset "Exponential" begin
                dist = Exponential()

                f = x -> exp(-x)
                ipdf = y -> -log(y)
                tailarea = x -> exp(-x)
                fallback_generator = x -> rng -> x - log1p(-rand(rng, T))
                domain = (T(0), T(Inf))

                @testset for N in [1, 2, 3, 4]
                    z = UnboundedZiggurat(f, domain, N; ipdf, tailarea, fallback_generator)
                    @test typeof(rand(z)) == eltype(z) == ZigguratTools.Ytype(z) == T
                    test_samples(z, dist)
                end
            end

            @testset "SteppedExponential" begin
                dist = SteppedExponential()

                f = Base.Fix1(pdf, dist)
                ta = Base.Fix1(ccdf, dist)
                domain = (T(0), T(Inf))

                @testset for N in [1, 2, 3, 4]
                    z = UnboundedZiggurat(f, domain, N; tailarea = ta)
                    @test typeof(rand(z)) == eltype(z) == ZigguratTools.Ytype(z) == T
                    test_samples(z, dist)
                end
            end
        end

        @testset "Bounded Ziggurats" begin
            @testset "Truncated Normal (0.5 <= x <= 1)" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> √(-2log(√T(2π)*y))
                z = BoundedZiggurat(f, (T(0.5), T(1)), 10; ipdf)

                @test typeof(rand(z)) == eltype(z) == Ytype(z) == T
                test_samples(z, truncated(dist; lower = 0.5, upper = 1))
            end

            @testset "Truncated Normal (-1 <= x <= -0.5)" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> -√(-2log(√T(2π)*y))
                z = BoundedZiggurat(f, (T(-1), T(-0.5)), 10; ipdf)

                @test typeof(rand(z)) == eltype(z) == Ytype(z) == T
                test_samples(z, truncated(dist; lower = -1, upper = -0.5))
            end

            @testset "Truncated Normal (0.5 <= x <= 10)" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> √(-2log(√T(2π)*y))
                z = BoundedZiggurat(f, (T(0.5), T(10)), 10; ipdf)

                @test typeof(rand(z)) == eltype(z) == Ytype(z) == T
                test_samples(z, truncated(dist; lower = 0.5, upper = 10))
            end

            @testset "Truncated Normal (-10 <= x <= -0.5)" begin
                dist = Normal()

                f = x -> 1/√T(2π) * exp(-x^2/2)
                ipdf = y -> -√(-2log(√T(2π)*y))
                z = BoundedZiggurat(f, (T(-10), T(-0.5)), 10; ipdf)

                @test typeof(rand(z)) == eltype(z) == Ytype(z) == T
                test_samples(z, truncated(dist; lower = -10, upper = -0.5))
            end
        end
    end
end

nothing
