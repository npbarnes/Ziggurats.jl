function test_common_layer_properties(x, y, N, modalboundary, slopesign, f, invf)
    # Initialization
    @test y[1] == zero(eltype(y))

    # Number of layers
    @test length(x) == length(y) == N + 1

    # Each layer has positive thickness
    @test all(y[i + 1] > y[i] for i in 1:N)

    # The ziggurat may extend to the left or right of the modal boundary, but not both.
    # The last x may be on the same side as the rest, or it may equal the modal boundary.
    lastsign = sign(modalboundary - x[end])
    @test lastsign == slopesign || lastsign == 0
    @test all(sign(modalboundary - x[i]) == slopesign for i in 1:N)

    # Ziggurats never get wider as you go up the tower (they can get narrower or
    # stay the same width)
    @test all(abs(x[i] - modalboundary) >= abs(x[i + 1] - modalboundary) for i in 1:N)

    # Layer areas are all equal
    A = abs(x[1] - modalboundary) * (y[2] - y[1])
    @test all(abs(x[i] - modalboundary) * (y[i + 1] - y[i]) ≈ A for i in 2:N)

    # An optimal ziggurat will have its upper boundary close to the pdf.
    @test y[end] ≈ pdf(modalboundary)

    # A valid ziggurat will have its upper boundary greater than or equal to the pdf.
    @test y[end] >= f(modalboundary)

    # x[i] is the generalized inverse of y[i] for all i in 2:N. i=1 is excluded
    # because distributions with an unbounded domain will have an artificial
    # base, and distributions with a bounded domain will always have
    # x[1]==argminboundary. i=N+1 is excluded because y might be greater
    # than pdf(mode), so ipdf(y) might not be defined.
    @test all(x[i] ≈ invf(y[i]) for i in 2:N)

    if y[end] == f(modalboundary)
        @test x[end] ≈ invf(y[end])
    end
end

function test_unbounded_domain(x, y, modalboundary, tailarea)
    true_base_area = abs(x[2] - modalboundary) * y[2] + tailarea(x[2])
    artificial_base_area = abs(x[1] - modalboundary) * y[2]

    @test true_base_area ≈ artificial_base_area
end

function test_bounded_domain(x, argminboundary)
    @test x[1] == argminboundary
end

function testset_body(
    N,
    modalboundary,
    argminboundary,
    slopesign,
    f,
    invf,
    tailarea;
    continuouspdf,
    initiallyflat
)
    if tailarea === nothing
        x, y = ZigguratTools.search(N, modalboundary, argminboundary, f, invf)
    else
        x, y = ZigguratTools.search(N, modalboundary, argminboundary, f, invf, tailarea)
    end

    # General tests
    test_common_layer_properties(x, y, N, modalboundary, slopesign, f, invf)

    # Additional tests for special cases
    if continuouspdf
        for i in 2:(N + 1)
            if x[i] != argminboundary
                @test y[i] ≈ f(x[i])
            end
        end
    end
    if !initiallyflat
        @test x[end] ≈ modalboundary atol = 1e-5
    end
    if tailarea === nothing
        test_bounded_domain(x, argminboundary)
    else
        test_unbounded_domain(x, y, modalboundary, tailarea)
    end

    x, y
end

# Test ziggurats produced from Distributions.jl.
function test_dist_ziggurats(
    Ns,
    dist,
    modalboundary,
    argminboundary;
    continuouspdf,
    initiallyflat,
    boundeddomain
)
    f = Base.Fix1(pdf, dist)
    domain = extrema(dist)
    invf = inverse(f, domain)

    @test domain == minmax(modalboundary, argminboundary)
    @test (modalboundary, argminboundary) == ZigguratTools._identify_mode(domain, f)

    if domain[1] == modalboundary # Decreasing distribution
        slopesign = -1
    elseif domain[2] == modalboundary # Increasing distribution
        slopesign = 1
    else
        error("modalboundary is not on the boundary of the support of dist.")
    end

    if boundeddomain
        tailarea = nothing
    else
        if slopesign == -1
            tailarea = Base.Fix1(ccdf, dist)
        elseif slopesign == 1
            tailarea = Base.Fix1(cdf, dist)
        else
            error("modalboundary is not on the boundary of the support of dist.")
        end
    end

    # Unnormalized
    uf(x) = 3 * f(x)
    uinvf(y) = invf(y / 3)
    if tailarea === nothing
        utailarea = nothing
    else
        utailarea = x -> 3 * tailarea(x)
    end

    @testset "Normalized pdf" begin
        @testset "N = $N" for N in Ns
            testset_body(
                N,
                modalboundary,
                argminboundary,
                slopesign,
                f,
                invf,
                tailarea;
                continuouspdf,
                initiallyflat
            )
        end
    end

    @testset "Unnormalized pdf" begin
        @testset "N = $N" for N in Ns
            testset_body(
                N,
                modalboundary,
                argminboundary,
                slopesign,
                uf,
                uinvf,
                utailarea;
                continuouspdf,
                initiallyflat
            )
        end
    end
end

@testset "Normal (x>=0)" begin
    dist = truncated(Normal(); lower = 0.0)
    test_dist_ziggurats(
        [1, 2, 256],
        dist,
        0.0,
        Inf;
        continuouspdf = true,
        initiallyflat = false,
        boundeddomain = false
    )
end

@testset "Normal (x<=0)" begin
    dist = truncated(Normal(); upper = 0.0)
    test_dist_ziggurats(
        [1, 2, 256],
        dist,
        0.0,
        -Inf;
        continuouspdf = true,
        initiallyflat = false,
        boundeddomain = false
    )
end

@testset "Non-standard Normal" begin
    dist = truncated(Normal(1.5, 3); upper = 1)
    test_dist_ziggurats(
        [1, 2, 256],
        dist,
        1.0,
        -Inf;
        continuouspdf = true,
        initiallyflat = false,
        boundeddomain = false
    )
end

@testset "Two Steps" begin
    # Designed so that the current build algorithm cannot produce an optimal ziggurat (i.e. y[end] ≈ f(mode))
    function f(x)
        if 0 <= x <= 1
            4.5
        elseif 1 < x <= 2
            1.0
        else
            0.0
        end
    end

    function invf(y)
        if 0 <= y <= 1
            2.0
        elseif 1 < y <= 4.5
            1.0
        else
            error()
        end
    end

    testset_body(
        3,
        0.0,
        2.0,
        -1,
        f,
        invf,
        nothing;
        continuouspdf = false,
        initiallyflat = true
    )
end

@testset "SteppedExponential" begin
    dist = SteppedExponential()
    test_dist_ziggurats(
        [1, 2, 256],
        dist,
        0.0,
        Inf;
        continuouspdf = false,
        initiallyflat = true,
        boundeddomain = false
    )
end

@testset "Truncated Normal (0.5 <= x <= 1)" begin
    dist = truncated(Normal(); lower = 0.5, upper = 1)
    test_dist_ziggurats(
        [1, 2, 256],
        dist,
        0.5,
        1.0;
        continuouspdf = true,
        initiallyflat = false,
        boundeddomain = true
    )
end

@testset "Truncated Normal (0.5 <= x <= 10)" begin
    dist = truncated(Normal(); lower = 0.5, upper = 10)
    test_dist_ziggurats(
        [1, 2, 256],
        dist,
        0.5,
        10.0;
        continuouspdf = true,
        initiallyflat = false,
        boundeddomain = true
    )
end

@testset "Truncated Normal (-1 <= x <= -0.5)" begin
    dist = truncated(Normal(); lower = -1, upper = -0.5)
    test_dist_ziggurats(
        [1, 2, 256],
        dist,
        -0.5,
        -1.0;
        continuouspdf = true,
        initiallyflat = false,
        boundeddomain = true
    )
end

@testset "Truncated Normal (-10 <= x <= -0.5)" begin
    dist = truncated(Normal(); lower = -10, upper = -0.5)
    test_dist_ziggurats(
        [1, 2, 256],
        dist,
        -0.5,
        -10.0;
        continuouspdf = true,
        initiallyflat = false,
        boundeddomain = true
    )
end

nothing
