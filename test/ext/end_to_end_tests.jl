@testset "End to End" begin
    @testset "Monotonic" begin
        @testset "Bounded" begin
            @testset "Upward Slope" begin
                dists = [
                    truncated(Normal(); lower = -3, upper = 0),
                    truncated(Normal(3, 5); lower = -1, upper = 1.5),
                    truncated(Cosine(0, 1); lower = 0)
                ]

                @testset for d in dists
                    z = monotonic_ziggurat(d)
                    @test z isa BoundedZiggurat
                    test_samples(z, d)
                end
            end

            @testset "Downward Slope" begin
                dists = [
                    truncated(Normal(); lower = 0, upper = 3),
                    truncated(Exponential(); upper = 3),
                    LogUniform(1, 4),
                    truncated(Semicircle(1); lower = 0)
                ]

                @testset for d in dists
                    z = monotonic_ziggurat(d)
                    @test z isa BoundedZiggurat
                    test_samples(z, d)
                end
            end
        end

        @testset "Unbounded" begin
            @testset "Upward Slope" begin
                dists =
                    [truncated(Normal(); upper = 0.0), truncated(Normal(3, 5); upper = 1.5)]

                @testset for d in dists
                    z = monotonic_ziggurat(d)
                    @test z isa UnboundedZiggurat
                    test_samples(z, d)
                end
            end

            @testset "Downward Slope" begin
                dists = [
                    truncated(Normal(); lower = 0.0),
                    truncated(Normal(-3, 5); lower = -1.5),
                    Exponential(),
                    SteppedExponential()
                ]

                @testset for d in dists
                    z = monotonic_ziggurat(d)

                    @test z isa UnboundedZiggurat
                    test_samples(z, d)
                end
            end
        end
    end

    @testset "Symmetric" begin
        @testset "Bounded" begin end

        @testset "Unbounded" begin end
    end

    @testset "Composite" begin
        @testset "Bounded" begin end

        @testset "Bounded Below" begin end

        @testset "Bounded Above" begin end

        @testset "Unbounded" begin end
    end
end
