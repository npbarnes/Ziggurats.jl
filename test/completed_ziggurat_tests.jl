function test_common_layer_properties(dist, x, y, N, modalboundary, tailarea, mypdf, myipdf)
    true_base_area = abs(x[2] - modalboundary) * y[2] + tailarea(x[2])
    artificial_base_area = abs(x[1] - modalboundary) * y[2]

    # Base area
    @test true_base_area ≈ artificial_base_area

    # boundary is mode
    @test modalboundary ≈ only(modes(dist))

    # Initialization
    @test y[1] == zero(eltype(y))

    # Number of layers
    @test length(x) == length(y) == N + 1

    # Each layer has positive thickness
    @test all(y[i + 1] > y[i] for i in 1:N)

    # The ziggurat may extend to the left or right of the modal boundary, but not both.
    # The last x may be on the same side as the rest, or it may equal the modal boundary.
    firstsign = sign(x[1] - modalboundary)
    lastsign = sign(x[N + 1] - modalboundary)
    @testset begin
        @test lastsign == firstsign || lastsign == 0
        @test all(sign(x[i] - modalboundary) == firstsign for i in 2:N)
    end

    # Ziggurats never get wider as you go up the tower (they can get narrower or
    # stay the same width)
    @test all(abs(x[i] - modalboundary) >= abs(x[i + 1] - modalboundary) for i in 1:N)

    # Layer areas are all equal
    @test all(abs(x[i] - modalboundary) * (y[i + 1] - y[i]) ≈ true_base_area for i in 1:N)

    # Upper boundary is close to and greater than or equal to f(m)
    @test y[end] ≈ mypdf(mode(dist)) && y[end] >= mypdf(mode(dist))

    # x = f^-1(y), but not for i=1 because of the artificial base and not for
    # i=N+1 because y might be greater than pdf(mode), so ipdf(y) might not be
    # defined.
    @test all(x[i] ≈ myipdf(y[i]) for i in 2:N)
end

function testset_body(
    N,
    dist,
    modalboundary,
    tailarea,
    f,
    invf;
    continuouspdf,
    initiallyflat
)
    x, y = ZigguratTools.search(N, modalboundary, tailarea, f, invf)
    test_common_layer_properties(dist, x, y, N, modalboundary, tailarea, f, invf)

    if continuouspdf
        @test all(y[i] ≈ f(x[i]) for i in 2:(N + 1))
    end
    if !initiallyflat
        @test x[end] ≈ modalboundary
    end

    x, y
end

function test_dist_ziggurats(Ns, dist, modalboundary; continuouspdf, initiallyflat)
    f = Base.Fix1(pdf, dist)

    if f(mode(dist)) != f(modalboundary)
        error("Incorrect modalboundary.")
    end

    L, R = minimum(dist), maximum(dist)

    if L == modalboundary
        # Right handed ziggurat
        invf = Base.Fix1(ipdf_right, dist)
        tailarea = Base.Fix1(ccdf, dist)
    elseif R == modalboundary
        # Left handed ziggurat
        invf = Base.Fix1(ipdf_left, dist)
        tailarea = Base.Fix1(cdf, dist)
    else
        error("modalboundary is not on the boundary of the support of dist.")
    end

    # Unnormalized
    uf(x) = 3 * f(x)
    uinvf(y) = invf(y / 3)
    utailarea(x) = 3 * tailarea(x)

    @testset "Normalized pdf" begin
        @testset "N = $N" for N in Ns
            testset_body(
                N,
                dist,
                modalboundary,
                tailarea,
                f,
                invf;
                continuouspdf,
                initiallyflat
            )
        end
    end

    @testset "Unnormalized pdf" begin
        @testset "N = $N" for N in Ns
            testset_body(
                N,
                dist,
                modalboundary,
                utailarea,
                uf,
                uinvf;
                continuouspdf,
                initiallyflat
            )
        end
    end
end

@testset "Normal (x>=0)" begin
    dist = truncated(Normal(); lower = 0.0)
    modalboundary = 0.0
    test_dist_ziggurats(
        [1, 256],
        dist,
        modalboundary;
        continuouspdf = true,
        initiallyflat = false
    )
end

@testset "Normal (x<=0)" begin
    dist = truncated(Normal(); upper = 0.0)
    modalboundary = 0.0
    test_dist_ziggurats(
        [1, 256],
        dist,
        modalboundary;
        continuouspdf = true,
        initiallyflat = false
    )
end

@testset "Non-standard Normal" begin
    dist = truncated(Normal(1.5, 3); upper = 1)
    modalboundary = 1.0
    test_dist_ziggurats(
        [1, 256],
        dist,
        modalboundary;
        continuouspdf = true,
        initiallyflat = false
    )
end

@testset "SteppedExponential" begin
    dist = SteppedExponential()
    modalboundary = 0.0
    test_dist_ziggurats(
        [1, 256],
        dist,
        modalboundary;
        continuouspdf = false,
        initiallyflat = true
    )
end
