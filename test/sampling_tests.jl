function testsampling(dist, z)
    values = test_samples(z, dist, 10000)

    @test eltype(values) === eltype(dist)
    @test !any(isnan, values)
    @test !any(isinf, values)
end

@testset "Normal (x>=0)" begin
    dist = Normal()

    # Because of the choice of domain this ziggurat will actually be sampling
    # from truncated(Normal(), lower=0.0). A small number N is chosen so that
    # the fallback branch gets chosen with high probability and its confidence
    # intervals can be tested.
    N = 3
    f = let dist=dist; x -> pdf(dist, x) end
    ta = let dist=dist; x -> ccdf(dist,x) end
    domain = (0, Inf)
    z = UnboundedZiggurat(
        N,
        domain,
        f,
        inverse(f, domain),
        ta
    )

    testsampling(truncated(dist; lower = mode(dist)), z)
end

@testset "Normal (x<=0)" begin
    dist = Normal()

    # Because of the choice of a domain, this ziggurat will actually be sampling
    # from truncated(Normal(), upper=0.0). A small number N is chosen so that
    # the fallback branch gets chosen with high probability and its confidence
    # intervals can be tested.
    N = 3
    f = let dist=dist; x -> pdf(dist, x) end
    ta = let dist=dist; x -> cdf(dist,x) end
    domain = (-Inf, 0)
    z = UnboundedZiggurat(
        N,
        domain,
        f,
        inverse(f, domain),
        ta
    )

    testsampling(truncated(dist; upper = mode(dist)), z)
end

@testset "Exponential" begin
    dist = Exponential()

    N = 3
    f = let dist=dist; x -> pdf(dist, x) end
    ta = let dist=dist; x -> ccdf(dist,x) end
    domain = (0, Inf)
    z = UnboundedZiggurat(
        N,
        domain,
        f,
        inverse(f, domain),
        ta
    )

    testsampling(dist, z)
end

@testset "SteppedExponential" begin
    dist = SteppedExponential()

    N = 3
    f = let dist=dist; x -> pdf(dist, x) end
    ta = let dist=dist; x -> ccdf(dist,x) end
    domain = (0, Inf)
    z = UnboundedZiggurat(
        N,
        domain,
        f,
        inverse(f, domain),
        ta
    )

    testsampling(dist, z)
end
