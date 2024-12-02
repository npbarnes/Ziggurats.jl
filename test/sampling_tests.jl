function testsampling(dist, z)
    values = test_samples(z, dist, 10000)

    @test eltype(values) === eltype(dist)
    @test !any(isnan, values)
    @test !any(isinf, values)
end
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

            testsampling(truncated(dist; lower = mode(dist)), z)
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

            testsampling(truncated(dist; upper = mode(dist)), z)
        end

        @testset "Exponential" begin
            dist = Exponential()

            N = 3
            f = Base.Fix1(pdf, dist)
            ta = Base.Fix1(ccdf, dist)
            domain = (0, Inf)
            z = UnboundedZiggurat(f, domain, N; tailarea = ta)

            testsampling(dist, z)
        end

        @testset "SteppedExponential" begin
            dist = SteppedExponential()

            N = 3
            f = Base.Fix1(pdf, dist)
            ta = Base.Fix1(ccdf, dist)
            domain = (0, Inf)
            z = UnboundedZiggurat(f, domain, N; tailarea = ta)

            testsampling(dist, z)
        end
    end

    @testset "Bounded Ziggurats" begin
        @testset "Truncated Normal (0.5 <= x <= 1)" begin
            dist = Normal()
    
            f = Base.Fix1(pdf, dist)
            z = BoundedZiggurat(f, (0.5, 1), 10)
    
            testsampling(truncated(dist; lower=0.5, upper=1), z)
        end
    
        @testset "Truncated Normal (-1 <= x <= -0.5)" begin
            dist = Normal()
    
            f = Base.Fix1(pdf, dist)
            z = BoundedZiggurat(f, (-1, -0.5), 10)
    
            testsampling(truncated(dist; lower=-1, upper=-0.5), z)
        end
    
        @testset "Truncated Normal (0.5 <= x <= 10)" begin
            dist = Normal()
    
            f = Base.Fix1(pdf, dist)
            z = BoundedZiggurat(f, (0.5, 10), 10)
    
            testsampling(truncated(dist; lower=0.5, upper=10), z)
        end
    
        @testset "Truncated Normal (-10 <= x <= -0.5)" begin
            dist = Normal()
    
            f = Base.Fix1(pdf, dist)
            z = BoundedZiggurat(f, (-10, -0.5), 10)
    
            testsampling(truncated(dist; lower=-10, upper=-0.5), z)
        end
    end
end

nothing
