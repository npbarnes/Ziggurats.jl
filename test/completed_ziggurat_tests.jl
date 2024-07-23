
using ZigguratTools, Test
using Distributions

function test_layer_properties(dist, N, z::UnboundedDecreasingZiggurat)
    x = z.x
    y = z.y
    A = z.layerarea

    # Number of layers
    @test length(x) == length(y) == N

    # Each layer has positive thickness
    @test y[begin] > zero(y[begin])
    @test all(y[i+1] > y[i] for i in eachindex(y)[begin:end-1])

    # Ziggurats never get wider (they can get narrower or stay the same width)
    @test all(x[i+1] <= x[i] for i in eachindex(x)[begin:end-1])

    # y = f(x)
    @test all(y .≈ pdf.(dist, x))

    # Base layer area
    @test x[1]*y[1] + ccdf(dist, x[1]) ≈ A 

    # Other layer areas
    @test all(x[i] * (y[i+1] - y[i]) ≈ A for i in eachindex(x)[1:end-1]) 

    # Upper boundary is close to and greater than or equal to f(m)
    @test y[end] ≈ pdf(dist, mode(dist)) && y[end] >= pdf(dist, mode(dist))
end

@testset "Normal (x>=0)" begin
    dist = Normal()
    N = 256
    
    z = UnboundedDecreasingZiggurat(
        x->pdf(dist,x),
        y->ipdf_right(dist,y),
        x->ccdf(dist,x),
        mode(dist),
        N,
        x->sampler(truncated(dist, lower=x))
    )

    test_layer_properties(dist, N, z)
end

@testset "Exponential" begin
    dist = Exponential()
    N = 256

    z = UnboundedDecreasingZiggurat(
        x->pdf(dist,x),
        y->ipdf_right(dist,y),
        x->ccdf(dist,x),
        mode(dist),
        N,
        x->sampler(truncated(dist, lower=x))
    )

    test_layer_properties(dist, N, z)
end

@testset "SteppedExponential" begin
    include("./testutils.jl")

    dist = SteppedExponential()
    N = 256

    z = UnboundedDecreasingZiggurat(
        x->pdf(dist,x),
        y->ipdf_right(dist,y),
        x->ccdf(dist,x),
        mode(dist),
        N,
        x->sampler(truncated(dist, lower=x))
    )

    test_layer_properties(dist, N, z)
end