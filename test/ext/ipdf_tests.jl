@testset "Distribution ipdf" begin
    ipdf = DistributionsExt.ipdf
    @testset "y above pdf(mode) throws domain error" begin
        dists = [
            Exponential(),
            truncated(Normal(), lower=0),
            truncated(Normal(), upper=0)
        ]

        for d in dists
            @testset let d=d
                @test_throws DomainError ipdf(d, pdf(d, mode(d)) + 1)
            end
        end
    end

    @testset "y below 0 throws domain error" begin
        dists = [
            Exponential(),
            truncated(Normal(), lower=0),
            truncated(Normal(), upper=0)
        ]

        for d in dists
            @testset let d=d
                @test_throws DomainError ipdf(d, -1)
            end
        end
    end

    @testset "fail silently for most non-monotonic distributions" begin
        
    end
end
