@testset "JET.jl Tests" begin
    @testset for T in (Float64, Float32, Float16)
        @testset for N in [1, 2, 3, 4]
            @testset "BoundedZiggurats" begin
                @testset "Default args" begin
                    z = @inferred BoundedZiggurat{T,T,N}(x->exp(-x), (T(0), T(1)))
                    @test_opt rand(z)
                    @test_call rand(z)
                end
            end
            @testset "UnboundedZiggurats" begin
                @testset "Default args" begin
                    z = @inferred UnboundedZiggurat{T,T,N}(x->exp(-x), (T(0), typemax(T)))
                    @test_opt rand(z)
                    @test_call rand(z)
                end
                @testset "Overriding tailarea" begin
                    z = @inferred UnboundedZiggurat{T,T,N}(x->exp(-x), (T(0), typemax(T)); tailarea = x->exp(-x))
                    @test_opt rand(z)
                    @test_call rand(z)
                end
                @testset "Overriding fallback" begin
                    z = @inferred UnboundedZiggurat{T,T,N}(
                        x->exp(-x),
                        (T(0), typemax(T));
                        fallback = (rng, a) -> a-log1p(-rand(rng, typeof(a)))
                    )
                    @test_opt rand(z)
                    @test_call rand(z)
                end
            end
        end
    end
end
