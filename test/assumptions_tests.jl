@testset "Misc assumptions" begin
    @testset "pdf is never evaluated at infinity" begin
        f = x -> isinf(x) ? error("pdf evaluated at infinity") : exp(-x)

        @testset "PDFWrap doesn't evaluate pdf at infinity" begin
            wf = Ziggurats.PDFWrap(f, 0, Inf)
            @test wf(Inf) == 0
        end

        @testset "end to end, pdf is not evaluated at infinity" begin
            @testset for T in TestTypes, N in [1, 2, 3, 4]
                z = ziggurat(f, (T(0), T(Inf)), N)
                @test rand(z, 1000) isa Any
            end
        end
    end

    @testset "rand(z) is type stable even when pdf is not" begin
        f = x -> x < 0.5 ? 1 : 0.5
        @testset for T in TestTypes, N in [1, 2, 3, 4]
            z = ziggurat(f, (T(0), T(1)), N)
            # Type stable, not type grounded, so we use @inferred instead of @report_opt
            @test @inferred(rand(z)) isa T
        end
    end

    @testset "integers get promoted to floats wherever necessary" begin
        f = x -> x < 1 ? 2 : 1 # evaluates to Int
        invf = y -> 1 <= y <= 2 ? 1 : 2 # evaluates to Int
        @testset for N in [1, 2, 3, 4], ipdf in [nothing, invf]
            z = ziggurat(f, (0, 2), N; ipdf)
            @test eltype(z) <: AbstractFloat
            @test Ziggurats.Ytype(z) <: AbstractFloat
            @test @inferred(rand(z)) isa AbstractFloat
        end
    end
end

nothing
