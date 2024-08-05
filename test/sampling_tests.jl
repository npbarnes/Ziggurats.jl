function testsampling(dist, z)
    values = test_samples(z, dist, 10000)

    @test eltype(values) === eltype(dist)
    @test !any(isnan, values)
    @test !any(isinf, values)
end

@testset "Normal (x>=0)" begin
    dist = Normal()

    # Because of the choice of UnboundedDecreasingZiggurat, ipdf_right, and ccdf,
    # this ziggurat will actually be sampling from truncated(Normal(), lower=0.0).
    # A small number N is chosen so that the fallback branch gets chosen with high
    # probability and its confidence intervals can be tested.
    N = 3
    z = monotonic_ziggurat(
        N,
        mode(dist),
        x -> ccdf(dist, x),
        x -> pdf(dist, x),
        y -> ipdf_right(dist, y),
        x -> sampler(truncated(dist; lower = x))
    )

    testsampling(truncated(dist; lower = mode(dist)), z)
end

@testset "Normal (x<=0)" begin
    dist = Normal()

    # Because of the choice of a AbstractUnboundedMonotonicZiggurat, ipdf_left, and cdf,
    # this ziggurat will actually be sampling from truncated(Normal(), lower=0.0).
    # A small number N is chosen so that the fallback branch gets chosen with high
    # probability and its confidence intervals can be tested.
    N = 3
    z = monotonic_ziggurat(
        N,
        mode(dist),
        x -> cdf(dist, x),
        x -> pdf(dist, x),
        y -> ipdf_left(dist, y),
        x -> sampler(truncated(dist; upper = x))
    )

    testsampling(truncated(dist; upper = mode(dist)), z)
end

@testset "Exponential" begin
    dist = Exponential()

    N = 3
    z = monotonic_ziggurat(
        N,
        mode(dist),
        x -> ccdf(dist, x),
        x -> pdf(dist, x),
        y -> ipdf_right(dist, y),
        x -> sampler(truncated(dist; lower = x))
    )

    testsampling(dist, z)
end

@testset "SteppedExponential" begin
    dist = SteppedExponential()

    N = 3
    z = monotonic_ziggurat(
        N,
        mode(dist),
        x -> ccdf(dist, x),
        x -> pdf(dist, x),
        y -> ipdf_right(dist, y),
        x -> sampler(truncated(dist; lower = x))
    )

    testsampling(dist, z)
end
