@testset "Sampling Tests" begin
    @testset "Unbounded Ziggurats" begin
        @testset "Normal (x>=0)" begin
            dist = Normal()

            # Because of the choice of domain this ziggurat will actually be sampling
            # from truncated(Normal(), lower=0.0). A small number N is chosen so that
            # the fallback branch gets chosen with high probability and its confidence
            # intervals can be tested.
            N = 3
            f = Base.Fix1(pdf, dist)
            ta = Base.Fix1(ccdf, dist)
            domain = (0, Inf)
            z = UnboundedZiggurat(f, domain, N; tailarea = ta)

            test_samples(z, truncated(dist; lower = mode(dist)))
        end

        @testset "Normal (x<=0)" begin
            dist = Normal()

            # Because of the choice of a domain, this ziggurat will actually be sampling
            # from truncated(Normal(), upper=0.0). A small number N is chosen so that
            # the fallback branch gets chosen with high probability and its confidence
            # intervals can be tested.
            N = 3
            f = Base.Fix1(pdf, dist)
            ta = Base.Fix1(cdf, dist)
            domain = (-Inf, 0)
            z = UnboundedZiggurat(f, domain, N; tailarea = ta)

            test_samples(z, truncated(dist; upper = mode(dist)))
        end

        @testset "Exponential" begin
            dist = Exponential()

            N = 3
            f = Base.Fix1(pdf, dist)
            ta = Base.Fix1(ccdf, dist)
            domain = (0, Inf)
            z = UnboundedZiggurat(f, domain, N; tailarea = ta)

            test_samples(z, dist)
        end

        @testset "SteppedExponential" begin
            dist = SteppedExponential()

            N = 3
            f = Base.Fix1(pdf, dist)
            ta = Base.Fix1(ccdf, dist)
            domain = (0, Inf)
            z = UnboundedZiggurat(f, domain, N; tailarea = ta)

            test_samples(z, dist)
        end
    end

    @testset "Bounded Ziggurats" begin
        @testset "Truncated Normal (0.5 <= x <= 1)" begin
            dist = Normal()

            f = Base.Fix1(pdf, dist)
            z = BoundedZiggurat(f, (0.5, 1), 10)

            test_samples(z, truncated(dist; lower = 0.5, upper = 1))
        end

        @testset "Truncated Normal (-1 <= x <= -0.5)" begin
            dist = Normal()

            f = Base.Fix1(pdf, dist)
            z = BoundedZiggurat(f, (-1, -0.5), 10)

            test_samples(z, truncated(dist; lower = -1, upper = -0.5))
        end

        @testset "Truncated Normal (0.5 <= x <= 10)" begin
            dist = Normal()

            f = Base.Fix1(pdf, dist)
            z = BoundedZiggurat(f, (0.5, 10), 10)

            test_samples(z, truncated(dist; lower = 0.5, upper = 10))
        end

        @testset "Truncated Normal (-10 <= x <= -0.5)" begin
            dist = Normal()

            f = Base.Fix1(pdf, dist)
            z = BoundedZiggurat(f, (-10, -0.5), 10)

            test_samples(z, truncated(dist; lower = -10, upper = -0.5))
        end
    end
end

nothing
