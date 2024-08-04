function test_continuous_distribution_layers(
    dist,
    x,
    y,
    N,
    modalboundary,
    tailarea,
    mypdf,
    myipdf
)
    # y = f(x), except for the artificial base
    @test all(y[i] ≈ mypdf(x[i]) for i in 2:(N + 1))
end

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
    @test lastsign == firstsign || lastsign == 0
    @test all(sign(x[i] - modalboundary) == firstsign for i in 2:N)

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

@testset "Normal (x>=0)" begin
    dist = truncated(Normal(); lower = 0.0)
    modalboundary = 0.0
    @testset "Normalized pdf" begin
        tailarea = x -> ccdf(dist, x)
        mypdf = x -> pdf(dist, x)
        myipdf = y -> ipdf_right(dist, y)
        @testset "N = 256" begin
            N = 256
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
            test_continuous_distribution_layers(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
        @testset "N = 1" begin
            N = 1
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
            test_continuous_distribution_layers(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
    end
    @testset "Unnormalized pdf" begin
        tailarea = x -> ccdf(dist, x)
        mypdf = x -> 3 * pdf(dist, x)
        myipdf = y -> ipdf_right(dist, y / 3)
        @testset "N = 256" begin
            N = 256
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
            test_continuous_distribution_layers(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
        @testset "N = 1" begin
            N = 1
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
            test_continuous_distribution_layers(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
    end
end

@testset "Normal (x<=0)" begin
    dist = truncated(Normal(); upper = 0.0)
    modalboundary = 0.0
    @testset "Normalized pdf" begin
        tailarea = x -> cdf(dist, x)
        mypdf = x -> pdf(dist, x)
        myipdf = y -> ipdf_left(dist, y)
        @testset "N = 256" begin
            N = 256
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
            test_continuous_distribution_layers(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
        @testset "N = 1" begin
            N = 1
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
            test_continuous_distribution_layers(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
    end
    @testset "Unnormalized pdf" begin
        tailarea = x -> cdf(dist, x)
        mypdf = x -> 3 * pdf(dist, x)
        myipdf = y -> ipdf_left(dist, y / 3)
        @testset "N = 256" begin
            N = 256
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
            test_continuous_distribution_layers(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
        @testset "N = 1" begin
            N = 1
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
            test_continuous_distribution_layers(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
    end
end

@testset "SteppedExponential" begin
    dist = SteppedExponential()
    modalboundary = 0.0
    @testset "Normalized pdf" begin
        tailarea = x -> ccdf(dist, x)
        mypdf = x -> pdf(dist, x)
        myipdf = y -> ipdf_right(dist, y)
        @testset "N = 256" begin
            N = 256
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
        @testset "N = 1" begin
            N = 1
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
    end
    @testset "Unnormalized pdf" begin
        tailarea = x -> ccdf(dist, x)
        mypdf = x -> 3 * pdf(dist, x)
        myipdf = y -> ipdf_right(dist, y / 3)
        @testset "N = 256" begin
            N = 256
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
        @testset "N = 1" begin
            N = 1
            x, y = ZigguratTools.search(N, modalboundary, tailarea, mypdf, myipdf)

            test_common_layer_properties(
                dist,
                x,
                y,
                N,
                modalboundary,
                tailarea,
                mypdf,
                myipdf
            )
        end
    end
end