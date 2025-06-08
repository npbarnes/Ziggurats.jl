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

            @test z_call.x == z_lamb.x == z_func.x
            @test z_call.y == z_lamb.y == z_func.y
            @test z_call.modalboundary == z_lamb.modalboundary == z_func.modalboundary
        end

        @testset "UnboundedZiggurat" begin
            z_call = UnboundedZiggurat(call, (0, Inf), 8)
            z_lamb = UnboundedZiggurat(lamb, (0, Inf), 8)
            z_func = UnboundedZiggurat(func, (0, Inf), 8)

            @test z_call.x == z_lamb.x == z_func.x
            @test z_call.y == z_lamb.y == z_func.y
            @test z_call.modalboundary == z_lamb.modalboundary == z_func.modalboundary
        end
    end

    @testset "Argument checking utilities" begin
        @testset "_choose_tailarea_func" begin
            inc_func(x) = x
            dec_func(x) = -x

            notnothing1 = x->1
            notnothing2 = x->2
            notnothing3 = x->3
            
            #! format: off
            @test                               ZigguratTools._choose_tailarea_func(inc_func, [0,1], nothing, nothing, nothing) === nothing
            @test                               ZigguratTools._choose_tailarea_func(dec_func, [0,1], nothing, nothing, nothing) === nothing

            @test_throws "a ccdf is provided"   ZigguratTools._choose_tailarea_func(inc_func, [0,1], nothing, nothing, notnothing1)
            @test                               ZigguratTools._choose_tailarea_func(dec_func, [0,1], nothing, nothing, notnothing1) === notnothing1

            @test                               ZigguratTools._choose_tailarea_func(inc_func, [0,1], nothing, notnothing1, nothing) === notnothing1
            @test_throws "a cdf is provided"    ZigguratTools._choose_tailarea_func(dec_func, [0,1], nothing, notnothing1, nothing)

            @test                               ZigguratTools._choose_tailarea_func(inc_func, [0,1], nothing, notnothing1, notnothing2) === notnothing1
            @test                               ZigguratTools._choose_tailarea_func(dec_func, [0,1], nothing, notnothing1, notnothing2) === notnothing2

            @test                               ZigguratTools._choose_tailarea_func(inc_func, [0,1], notnothing1, nothing, nothing) === notnothing1
            @test                               ZigguratTools._choose_tailarea_func(dec_func, [0,1], notnothing1, nothing, nothing) === notnothing1

            @test                               ZigguratTools._choose_tailarea_func(inc_func, [0,1], notnothing1, nothing, notnothing2) === notnothing1
            @test                               ZigguratTools._choose_tailarea_func(dec_func, [0,1], notnothing1, nothing, notnothing2) === notnothing1
            
            @test                               ZigguratTools._choose_tailarea_func(inc_func, [0,1], notnothing1, notnothing2, nothing) === notnothing1
            @test                               ZigguratTools._choose_tailarea_func(dec_func, [0,1], notnothing1, notnothing2, nothing) === notnothing1

            @test                               ZigguratTools._choose_tailarea_func(inc_func, [0,1], notnothing1, notnothing2, notnothing3) === notnothing1
            @test                               ZigguratTools._choose_tailarea_func(dec_func, [0,1], notnothing1, notnothing2, notnothing3) === notnothing1
            #! format: on
        end
    end
end
