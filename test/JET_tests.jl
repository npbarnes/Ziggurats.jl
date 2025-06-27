@testset "JET.jl Tests" begin
    @testset for T in (Float64, Float32, Float16)
        @testset for N in [1, 2, 3, 4]
            @testset "BoundedZiggurats" begin
                @testset "Default args" begin
                    z = BoundedZiggurat(x->exp(-x), (T(0), T(1)), N)
                    @test_opt rand(z)
                    @test_call rand(z)
                end
            end
            @testset "UnboundedZiggurats" begin
                # TODO: Instead of target_modules=(Ziggurats,), use a filter that prunes
                # only specific parts of the call tree that are known to be problematic.
                # I think all the errors detected by JET ultimately come from the `verbose`
                # flag in Roots.jl, but in normal use verbose is always false. It's only for
                # debugging, so I don't care if there's runtime dispatch on that branch.
                # I want to exclude that particular line of code in Roots.jl, but I don't
                # know how to get JET to do that.
                @testset "Default args" begin
                    z = UnboundedZiggurat(x->exp(-x), (T(0), typemax(T)), N)
                    @test_opt target_modules=(Ziggurats,) rand(z)
                    @test_call rand(z)
                end
                @testset "Overriding tailarea" begin
                    z = UnboundedZiggurat(
                        x->exp(-x),
                        (T(0), typemax(T)),
                        N;
                        tailarea = x->exp(-x)
                    )
                    @test_opt target_modules=(Ziggurats,) rand(z)
                    @test_call rand(z)
                end
                @testset "Overriding fallback" begin
                    z = UnboundedZiggurat(
                        x->exp(-x),
                        (T(0), typemax(T)),
                        N;
                        fallback = (rng, a) -> a-log1p(-rand(rng, typeof(a)))
                    )
                    @test_opt target_modules=(Ziggurats,) rand(z)
                    @test_call rand(z)
                end
            end
        end
    end
end
