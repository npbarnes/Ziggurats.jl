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

    A = abs(x[1] - modalboundary) * (y[2] - y[1])
    # Layer area is positive
    @test A > 0
    # Layer areas are all equal
    @test all(abs(x[i] - modalboundary) * (y[i + 1] - y[i]) ≈ A for i in 2:N)

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

function test_not_initially_flat(x, modalboundary; atol = 1e-8)
    @test x[end] ≈ modalboundary atol = atol
end

function test_continuous_pdf(x, y, f, N, argminboundary)
    # Expect f(x) = y, but not on a discontinuity. The argminboundary is a discontinuity for
    # bounded ziggurats, so it is excluded from this test.
    for i in 2:(N + 1)
        if x[i] != argminboundary
            @test y[i] ≈ f(x[i])
        end
    end
end

@testset "Completed Ziggurat" begin
    @testset "Domain type=$T" for T in (Float16, Float32, Float64)
        @testset "Normal (x>=0)" begin
            f = x -> 1/√T(2π) * exp(-x^2/2)
            invf = y -> √(-2log(√T(2π)*y))
            tailarea = x -> erfc(x/√2)/2

            uf(x) = 3 * f(x)
            uinvf(y) = √(-2log(√T(2π)*y/3))
            utailarea(x) = 3 * tailarea(x)

            modalboundary = T(0.0)
            argminboundary = T(Inf)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )
            wuf = ZigguratTools.PDFWrap(uf, modalboundary, argminboundary)
            wuinvf = ZigguratTools.IPDFWrap(
                uinvf,
                modalboundary,
                argminboundary,
                wuf(modalboundary),
                wuf(argminboundary)
            )

            slopesign = -1

            Ns = [1, 2, 3, 4]

            @testset "Normalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(
                        N,
                        modalboundary,
                        argminboundary,
                        wf,
                        winvf,
                        tailarea
                    )

                    @test eltype(x) == eltype(y) == T

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wf,
                        winvf
                    )

                    test_unbounded_domain(x, y, modalboundary, tailarea)
                    #test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary; atol = 1e-5)
                end
            end

            @testset "Unnormalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(
                        N,
                        modalboundary,
                        argminboundary,
                        wuf,
                        wuinvf,
                        utailarea
                    )

                    @test eltype(x) == eltype(y) == T

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wuf,
                        wuinvf
                    )

                    test_unbounded_domain(x, y, modalboundary, utailarea)
                    #test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wuf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary; atol = 1e-5)
                end
            end
        end

        @testset "Normal (x<=0)" begin
            f = x -> 1/√T(2π) * exp(-x^2/2)
            invf = y -> -√(-2log(√T(2π)*y))
            tailarea = x -> (1 + erf(x/√2))/2

            uf(x) = 3 * f(x)
            uinvf(y) = invf(y / 3)
            utailarea(x) = 3 * tailarea(x)

            modalboundary = T(0.0)
            argminboundary = T(-Inf)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )
            wuf = ZigguratTools.PDFWrap(uf, modalboundary, argminboundary)
            wuinvf = ZigguratTools.IPDFWrap(
                uinvf,
                modalboundary,
                argminboundary,
                wuf(modalboundary),
                wuf(argminboundary)
            )

            slopesign = 1

            Ns = [1, 2, 3, 4]

            @testset "Normalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(
                        N,
                        modalboundary,
                        argminboundary,
                        wf,
                        winvf,
                        tailarea
                    )

                    @test eltype(x) == eltype(y) == T

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wf,
                        winvf
                    )

                    test_unbounded_domain(x, y, modalboundary, tailarea)
                    #test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary; atol = 1e-5)
                end
            end

            @testset "Unnormalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(
                        N,
                        modalboundary,
                        argminboundary,
                        wuf,
                        wuinvf,
                        utailarea
                    )

                    @test eltype(x) == eltype(y) == T

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wuf,
                        wuinvf
                    )

                    test_unbounded_domain(x, y, modalboundary, utailarea)
                    #test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wuf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary; atol = 1e-5)
                end
            end
        end

        @testset "Non-standard Normal" begin
            μ = T(1.5)
            σ = T(3)
            f = x -> 1/(σ*√T(2π)) * exp(-(x-μ)^2/2σ)
            invf = y -> -√(-2σ * log(σ*√T(2π)*y)) + μ
            tailarea = x -> (1 + erf((x-μ)/(σ*√2)))/2

            uf(x) = 3 * f(x)
            uinvf(y) = invf(y / 3)
            utailarea(x) = 3 * tailarea(x)

            modalboundary = T(1.0)
            argminboundary = T(-Inf)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )
            wuf = ZigguratTools.PDFWrap(uf, modalboundary, argminboundary)
            wuinvf = ZigguratTools.IPDFWrap(
                uinvf,
                modalboundary,
                argminboundary,
                wuf(modalboundary),
                wuf(argminboundary)
            )

            slopesign = 1

            Ns = [1, 2, 3, 4]

            @testset "Normalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(
                        N,
                        modalboundary,
                        argminboundary,
                        wf,
                        winvf,
                        tailarea
                    )

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wf,
                        winvf
                    )

                    test_unbounded_domain(x, y, modalboundary, tailarea)
                    #test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end

            @testset "Unnormalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(
                        N,
                        modalboundary,
                        argminboundary,
                        wuf,
                        wuinvf,
                        utailarea
                    )

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wuf,
                        wuinvf
                    )

                    test_unbounded_domain(x, y, modalboundary, utailarea)
                    #test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wuf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end
        end

        @testset "Two Steps" begin
            # Designed so that the current build algorithm cannot produce an optimal ziggurat (i.e. y[end] ≈ f(mode))
            function f(x)
                if 0 <= x <= 1
                    T(4.5)
                elseif 1 < x <= 2
                    T(1.0)
                else
                    T(0.0)
                end
            end

            function invf(y)
                if 0 <= y <= 1
                    T(2.0)
                elseif 1 < y <= 4.5
                    T(1.0)
                else
                    error()
                end
            end

            N = 3
            modalboundary = T(0.0)
            argminboundary = T(2.0)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )

            slopesign = -1

            x, y = ZigguratTools.search(N, modalboundary, argminboundary, wf, winvf)

            test_common_layer_properties(x, y, N, modalboundary, slopesign, wf, winvf)

            #test_unbounded_domain(x, y, modalboundary, tailarea)
            test_bounded_domain(x, argminboundary)

            #test_continuous_pdf(x, y, f, N, argminboundary)
            #test_not_initially_flat(x, modalboundary)
        end

        @testset "SteppedExponential" begin
            dist = SteppedExponential(oneunit(T))
            f = Base.Fix1(pdf, dist)
            invf = inverse(f, extrema(dist)) # TODO: replace with direct implementation of inverse
            tailarea = Base.Fix1(ccdf, dist)

            uf(x) = 3 * f(x)
            uinvf(y) = invf(y / 3)
            utailarea(x) = 3 * tailarea(x)

            modalboundary = T(0.0)
            argminboundary = T(Inf)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )
            wuf = ZigguratTools.PDFWrap(uf, modalboundary, argminboundary)
            wuinvf = ZigguratTools.IPDFWrap(
                uinvf,
                modalboundary,
                argminboundary,
                wuf(modalboundary),
                wuf(argminboundary)
            )

            slopesign = -1

            Ns = [1, 2, 3, 4]

            @testset "Normalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(
                        N,
                        modalboundary,
                        argminboundary,
                        wf,
                        winvf,
                        tailarea
                    )

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wf,
                        winvf
                    )

                    test_unbounded_domain(x, y, modalboundary, tailarea)
                    #test_bounded_domain(x, argminboundary)

                    #test_continuous_pdf(x, y, f, N, argminboundary)
                    #test_not_initially_flat(x, modalboundary)
                end
            end

            @testset "Unnormalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(
                        N,
                        modalboundary,
                        argminboundary,
                        wuf,
                        wuinvf,
                        utailarea
                    )

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wuf,
                        wuinvf
                    )

                    test_unbounded_domain(x, y, modalboundary, utailarea)
                    #test_bounded_domain(x, argminboundary)

                    #test_continuous_pdf(x, y, uf, N, argminboundary)
                    #test_not_initially_flat(x, modalboundary)
                end
            end
        end

        @testset "Truncated Normal (0.5 <= x <= 1)" begin
            f = x -> 1/√T(2π) * exp(-x^2/2)
            invf = y -> √(-2log(√T(2π)*y))

            uf(x) = 3 * f(x)
            uinvf(y) = invf(y / 3)

            modalboundary = T(0.5)
            argminboundary = T(1.0)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )
            wuf = ZigguratTools.PDFWrap(uf, modalboundary, argminboundary)
            wuinvf = ZigguratTools.IPDFWrap(
                uinvf,
                modalboundary,
                argminboundary,
                wuf(modalboundary),
                wuf(argminboundary)
            )

            slopesign = -1

            Ns = [1, 2, 3, 4]

            @testset "Normalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(N, modalboundary, argminboundary, wf, winvf)

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wf,
                        winvf
                    )

                    #test_unbounded_domain(x, y, modalboundary, tailarea)
                    test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end

            @testset "Unnormalized pdf" begin
                @testset for N in Ns
                    x, y =
                        ZigguratTools.search(N, modalboundary, argminboundary, wuf, wuinvf)

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wuf,
                        wuinvf
                    )

                    #test_unbounded_domain(x, y, modalboundary, utailarea)
                    test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wuf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end
        end

        @testset "Truncated Normal (0.5 <= x <= 10)" begin
            f = x -> 1/√T(2π) * exp(-x^2/2)
            invf = y -> √(-2log(√T(2π)*y))

            uf(x) = 3 * f(x)
            uinvf(y) = invf(y / 3)

            modalboundary = T(0.5)
            argminboundary = T(10.0)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )
            wuf = ZigguratTools.PDFWrap(uf, modalboundary, argminboundary)
            wuinvf = ZigguratTools.IPDFWrap(
                uinvf,
                modalboundary,
                argminboundary,
                wuf(modalboundary),
                wuf(argminboundary)
            )

            slopesign = -1

            Ns = [1, 2, 3, 4]

            @testset "Normalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(N, modalboundary, argminboundary, wf, winvf)

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wf,
                        winvf
                    )

                    #test_unbounded_domain(x, y, modalboundary, tailarea)
                    test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end

            @testset "Unnormalized pdf" begin
                @testset for N in Ns
                    x, y =
                        ZigguratTools.search(N, modalboundary, argminboundary, wuf, wuinvf)

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wuf,
                        wuinvf
                    )

                    #test_unbounded_domain(x, y, modalboundary, utailarea)
                    test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wuf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end
        end

        @testset "Truncated Normal (-1 <= x <= -0.5)" begin
            f = x -> 1/√T(2π) * exp(-x^2/2)
            invf = y -> -√(-2log(√T(2π)*y))

            uf(x) = 3 * f(x)
            uinvf(y) = invf(y / 3)

            modalboundary = T(-0.5)
            argminboundary = T(-1.0)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )
            wuf = ZigguratTools.PDFWrap(uf, modalboundary, argminboundary)
            wuinvf = ZigguratTools.IPDFWrap(
                uinvf,
                modalboundary,
                argminboundary,
                wuf(modalboundary),
                wuf(argminboundary)
            )

            slopesign = 1

            Ns = [1, 2, 3, 4]

            @testset "Normalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(N, modalboundary, argminboundary, wf, winvf)

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wf,
                        winvf
                    )

                    #test_unbounded_domain(x, y, modalboundary, tailarea)
                    test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end

            @testset "Unnormalized pdf" begin
                @testset for N in Ns
                    x, y =
                        ZigguratTools.search(N, modalboundary, argminboundary, wuf, wuinvf)

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wuf,
                        wuinvf
                    )

                    #test_unbounded_domain(x, y, modalboundary, utailarea)
                    test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wuf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end
        end

        @testset "Truncated Normal (-10 <= x <= -0.5)" begin
            f = x -> 1/√T(2π) * exp(-x^2/2)
            invf = y -> -√(-2log(√T(2π)*y))

            uf(x) = 3 * f(x)
            uinvf(y) = invf(y / 3)

            modalboundary = T(-0.5)
            argminboundary = T(-10.0)

            wf = ZigguratTools.PDFWrap(f, modalboundary, argminboundary)
            winvf = ZigguratTools.IPDFWrap(
                invf,
                modalboundary,
                argminboundary,
                wf(modalboundary),
                wf(argminboundary)
            )
            wuf = ZigguratTools.PDFWrap(uf, modalboundary, argminboundary)
            wuinvf = ZigguratTools.IPDFWrap(
                uinvf,
                modalboundary,
                argminboundary,
                wuf(modalboundary),
                wuf(argminboundary)
            )

            slopesign = 1

            Ns = [1, 2, 3, 4]

            @testset "Normalized pdf" begin
                @testset for N in Ns
                    x, y = ZigguratTools.search(N, modalboundary, argminboundary, wf, winvf)

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wf,
                        winvf
                    )

                    #test_unbounded_domain(x, y, modalboundary, tailarea)
                    test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end

            @testset "Unnormalized pdf" begin
                @testset for N in Ns
                    x, y =
                        ZigguratTools.search(N, modalboundary, argminboundary, wuf, wuinvf)

                    test_common_layer_properties(
                        x,
                        y,
                        N,
                        modalboundary,
                        slopesign,
                        wuf,
                        wuinvf
                    )

                    #test_unbounded_domain(x, y, modalboundary, utailarea)
                    test_bounded_domain(x, argminboundary)

                    test_continuous_pdf(x, y, wuf, N, argminboundary)
                    test_not_initially_flat(x, modalboundary)
                end
            end
        end
    end
end

nothing
