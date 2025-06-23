@testset "Interface Tests" begin
    @testset "Different types of callables are accepted" begin
        struct Callable end
        (c::Callable)(x) = exp(-x^2)

        call = Callable()
        lamb = x -> exp(-x^2)
        func(x) = exp(-x^2)

        @testset "BoundedZiggurat" begin
            z_call = BoundedZiggurat(call, (0, 1), 8)
            z_lamb = BoundedZiggurat(lamb, (0, 1), 8)
            z_func = BoundedZiggurat(func, (0, 1), 8)

            import ZigguratTools: widths, layerratios, heights, highside
            @test widths(z_call) == widths(z_lamb) == widths(z_func)
            @test layerratios(z_call) == layerratios(z_lamb) == layerratios(z_func)
            @test heights(z_call) == heights(z_lamb) == heights(z_func)
            @test highside(z_call) == highside(z_lamb) == highside(z_func)
        end

        @testset "UnboundedZiggurat" begin
            z_call = UnboundedZiggurat(call, (0, Inf), 8)
            z_lamb = UnboundedZiggurat(lamb, (0, Inf), 8)
            z_func = UnboundedZiggurat(func, (0, Inf), 8)

            import ZigguratTools: widths, layerratios, heights, highside
            @test widths(z_call) == widths(z_lamb) == widths(z_func)
            @test layerratios(z_call) == layerratios(z_lamb) == layerratios(z_func)
            @test heights(z_call) == heights(z_lamb) == heights(z_func)
            @test highside(z_call) == highside(z_lamb) == highside(z_func)
        end
    end

end
