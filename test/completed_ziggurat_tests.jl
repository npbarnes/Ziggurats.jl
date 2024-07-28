
using ZigguratTools, Test
using Distributions

slopesign(z::MonotonicZiggurat) = sign(z.x[2] - z.x[1])

function baselayerarea(dist, z::UnboundedMonotonicZiggurat)
    if slopesign(z) > 0
        tail = cdf(dist, z.x[2])
    else
        tail = ccdf(dist, z.x[2])
    end
    abs(z.x[2] - z.modalboundary) * z.y[2] + tail
end

function baselayerarea(dist, z::BoundedMonotonicZiggurat)
    abs(z.x[1] - z.modalboundary) * z.y[2]
end

function layerwidth(z, l::Integer)
    if l < 1 || l >= length(z.x)
        error("Valid layer number are integers from 1 to N.")
    end
    abs(z.x[l] - z.modalboundary)
end

function layerarea(z, l::Integer)
    if l < 1 || l >= length(z.x)
        error("Valid layer number are integers from 1 to N.")
    end
    abs(z.x[l] - z.modalboundary) * (z.y[l+1] - z.y[l]) 
end

function test_continuous_distribution_layers(dist, z)
    x = z.x
    y = z.y

    # y = f(x)
    @test all(y ≈ pdf(dist, x) for (y, x) in Iterators.drop(zip(y,x), 1))
end

function test_layer_properties(dist, N, z)
    x = z.x
    y = z.y
    A = baselayerarea(dist, z)

    # Initialization
    @test y[1] == zero(eltype(y))

    # Number of layers
    @test length(x) == length(y) == N + 1

    # Each layer has positive thickness
    @test all(y[i+1] > y[i] for i in 1:N)

    # Ziggurats never get wider (they can get narrower or stay the same width)
    @test all(layerwidth(z, l) >= layerwidth(z, l+1) for l in 1:N-1)

    # Layer areas are all equal This also tests that the artificial layer area
    # of the first layer (x[1]*y[2]) equals the true baselayerarea derived from
    # the tail area function
    @test all(layerarea(z, l) ≈ A for l in 1:N) 

    # Upper boundary is close to and greater than or equal to f(m)
    @test y[end] ≈ pdf(dist, mode(dist)) && y[end] >= pdf(dist, mode(dist))

    # x = f^-1(y)
    #@test all(x ≈ ipdf(y) for (x, y) in Iterators.drop(zip(x, y), 1))
end

@testset "Normal (x>=0)" begin
    dist = Normal()
    N = 256
    
    z = monotonic_ziggurat(
        N,
        mode(dist),
        x->ccdf(dist,x),
        x->pdf(dist,x),
        y->ipdf_right(dist,y),
        x->sampler(truncated(dist, lower=x))
    )

    test_layer_properties(dist, N, z)
    test_continuous_distribution_layers(dist, z)
end

@testset "Normal (x<=0)" begin
    dist = Normal()
    N = 256
    
    z = monotonic_ziggurat(
        N,
        mode(dist),
        x->cdf(dist,x),
        x->pdf(dist,x),
        y->ipdf_left(dist,y),
        x->sampler(truncated(dist, upper=x))
    )

    test_layer_properties(dist, N, z)
    test_continuous_distribution_layers(dist, z)
end

@testset "Exponential" begin
    dist = Exponential()
    N = 256

    z = monotonic_ziggurat(
        N,
        mode(dist),
        x->ccdf(dist,x),
        x->pdf(dist,x),
        y->ipdf_right(dist,y),
        x->sampler(truncated(dist, lower=x))
    )

    test_layer_properties(dist, N, z)
    test_continuous_distribution_layers(dist, z)
end

@testset "SteppedExponential" begin
    include("./testutils.jl")

    dist = SteppedExponential()
    N = 256

    z = monotonic_ziggurat(
        N,
        mode(dist),
        x->ccdf(dist,x),
        x->pdf(dist,x),
        y->ipdf_right(dist,y),
        x->sampler(truncated(dist, lower=x))
    )

    test_layer_properties(dist, N, z)
end