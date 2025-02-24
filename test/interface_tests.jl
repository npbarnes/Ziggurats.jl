@testset "Interface Tests" begin
    @testset "Different types of callables are accepted" begin
        struct Callable end
        (c::Callable)(x) = exp(-x^2)
        
        call = Callable()
        lamb = x -> exp(-x^2)
        func(x) = exp(-x^2)

        @testset "BoundedZiggurat" begin
            z_call = BoundedZiggurat(call, (0,1), 8)
            z_lamb = BoundedZiggurat(lamb, (0,1), 8)
            z_func = BoundedZiggurat(func, (0,1), 8)

            @test z_call.x == z_lamb.x == z_func.x
            @test z_call.y == z_lamb.y == z_func.y
            @test z_call.modalboundary == z_lamb.modalboundary == z_func.modalboundary
        end

        @testset "UnboundedZiggurat" begin
            z_call = UnboundedZiggurat(call, (0,Inf), 8)
            z_lamb = UnboundedZiggurat(lamb, (0,Inf), 8)
            z_func = UnboundedZiggurat(func, (0,Inf), 8)

            @test z_call.x == z_lamb.x == z_func.x
            @test z_call.y == z_lamb.y == z_func.y
            @test z_call.modalboundary == z_lamb.modalboundary == z_func.modalboundary
        end
    end
end
