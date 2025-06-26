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
                # Filter `display` out because Roots.jl uses it to display the Tracks object
                # when it's run in `verbose` mode which causes runtime dispatch. But verbose
                # is only used for debugging, so under normal use this runtime dispatch never
                # occurs.
                function_filter(@nospecialize f) = f !== Base.Multimedia.display;
                @testset "Default args" begin
                    z = UnboundedZiggurat(x->exp(-x), (T(0), typemax(T)), N)
                    @test_opt function_filter=function_filter rand(z)
                    @test_call rand(z)
                end
                @testset "Overriding tailarea" begin
                    z = UnboundedZiggurat(
                        x->exp(-x),
                        (T(0), typemax(T)),
                        N;
                        tailarea = x->exp(-x)
                    )
                    @test_opt function_filter=function_filter rand(z)
                    @test_call rand(z)
                end
                @testset "Overriding fallback" begin
                    z = UnboundedZiggurat(
                        x->exp(-x),
                        (T(0), typemax(T)),
                        N;
                        fallback = (rng, a) -> a-log1p(-rand(rng, typeof(a)))
                    )
                    @test_opt function_filter=function_filter rand(z)
                    @test_call rand(z)
                end
            end
        end
    end
end
