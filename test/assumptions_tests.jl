@testset "pdf is never evaluated at infinity" begin
    f = x -> isinf(x) ? error("pdf evaluated at infinity") : exp(-x)

    @testset "PDFWrap doesn't evaluate pdf at infinity" begin
        try
            wf = Ziggurats.PDFWrap(f, 0, Inf)
        catch
            result = false
        else
            result = true
        end

        @test result

        try
            wf(Inf)
        catch
            result = false
        else
            result = true
        end

        @test result
    end

    @testset "end to end, pdf is not evaluated at infinity" begin
        @testset for T in TestTypes, N in [1, 2, 3, 4]
            try
                z = ziggurat(f, (T(0), T(Inf)), N)
                rand(z, 1000)
            catch
                result = false
            else
                result = true
            end

            @test result
        end
    end
end
