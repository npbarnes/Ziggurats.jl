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

    @testset "_check_arguments" begin
        @testset "Using a non-positive N throws an error" begin
            @testset "When domain is valid" begin
                @test_throws "N must" ZigguratTools._check_arguments(-1, (0, Inf))
                @test_throws "N must" ZigguratTools._check_arguments(-1, (0, 1))
                @test_throws "N must" ZigguratTools._check_arguments(-1, (-Inf, 1))
                @test_throws "N must" ZigguratTools._check_arguments(-1, (-Inf, Inf))

                @test_throws "N must" ZigguratTools._check_arguments(-1, [0, 1])
                @test_throws "N must" ZigguratTools._check_arguments(-1, [0, Inf])
                @test_throws "N must" ZigguratTools._check_arguments(-1, [-Inf, 1])
                @test_throws "N must" ZigguratTools._check_arguments(-1, [-Inf, Inf])

                @test_throws "N must" ZigguratTools._check_arguments(0, (0, 1))
                @test_throws "N must" ZigguratTools._check_arguments(0, (0, Inf))
                @test_throws "N must" ZigguratTools._check_arguments(0, (-Inf, 1))
                @test_throws "N must" ZigguratTools._check_arguments(0, (-Inf, Inf))

                @test_throws "N must" ZigguratTools._check_arguments(0, [0, 1])
                @test_throws "N must" ZigguratTools._check_arguments(0, [0, Inf])
                @test_throws "N must" ZigguratTools._check_arguments(0, [-Inf, 1])
                @test_throws "N must" ZigguratTools._check_arguments(0, [-Inf, Inf])

                @test_throws "N must" ZigguratTools._check_arguments(-1, (0, 0.5, 1))
                @test_throws "N must" ZigguratTools._check_arguments(-1, (0, 0.5, Inf))
                @test_throws "N must" ZigguratTools._check_arguments(-1, (-Inf, 0.5, 1))
                @test_throws "N must" ZigguratTools._check_arguments(-1, (-Inf, 0.5, Inf))

                @test_throws "N must" ZigguratTools._check_arguments(-1, [0, 0.5, 1])
                @test_throws "N must" ZigguratTools._check_arguments(-1, [0, 0.5, Inf])
                @test_throws "N must" ZigguratTools._check_arguments(-1, [-Inf, 0.5, 1])
                @test_throws "N must" ZigguratTools._check_arguments(-1, [-Inf, 0.5, Inf])

                @test_throws "N must" ZigguratTools._check_arguments(0, (0, 0.5, 1))
                @test_throws "N must" ZigguratTools._check_arguments(0, (0, 0.5, Inf))
                @test_throws "N must" ZigguratTools._check_arguments(0, (-Inf, 0.5, 1))
                @test_throws "N must" ZigguratTools._check_arguments(0, (-Inf, 0.5, Inf))

                @test_throws "N must" ZigguratTools._check_arguments(0, [0, 0.5, 1])
                @test_throws "N must" ZigguratTools._check_arguments(0, [0, 0.5, Inf])
                @test_throws "N must" ZigguratTools._check_arguments(0, [-Inf, 0.5, 1])
                @test_throws "N must" ZigguratTools._check_arguments(0, [-Inf, 0.5, Inf])
            end

            @testset "When the domain is invalid" begin
                @test_throws "N must be" ZigguratTools._check_arguments(-1, (0, 0))
                @test_throws "N must be" ZigguratTools._check_arguments(-1, (Inf, Inf))
                @test_throws "N must be" ZigguratTools._check_arguments(-1, (-Inf, -Inf))
                @test_throws "N must be" ZigguratTools._check_arguments(-1, (-Inf, Inf))

                @test_throws "N must be" ZigguratTools._check_arguments(0, (0, 0))
                @test_throws "N must be" ZigguratTools._check_arguments(0, (Inf, Inf))
                @test_throws "N must be" ZigguratTools._check_arguments(0, (-Inf, -Inf))
                @test_throws "N must be" ZigguratTools._check_arguments(0, (-Inf, Inf))

                @test_throws "N must be" ZigguratTools._check_arguments(-1, [0, 0])
                @test_throws "N must be" ZigguratTools._check_arguments(-1, [Inf, Inf])
                @test_throws "N must be" ZigguratTools._check_arguments(-1, [-Inf, -Inf])
                @test_throws "N must be" ZigguratTools._check_arguments(-1, [-Inf, Inf])

                @test_throws "N must be" ZigguratTools._check_arguments(0, [0, 0])
                @test_throws "N must be" ZigguratTools._check_arguments(0, [Inf, Inf])
                @test_throws "N must be" ZigguratTools._check_arguments(0, [-Inf, -Inf])
                @test_throws "N must be" ZigguratTools._check_arguments(0, [-Inf, Inf])
            end
        end

        # Empty domian
        @testset "Empty domains throw an empty domain error" begin
            @test_throws "empty domains" ZigguratTools._check_arguments(256, (0, 0))
            @test_throws "empty domains" ZigguratTools._check_arguments(256, (Inf, Inf))
            @test_throws "empty domains" ZigguratTools._check_arguments(256, (-Inf, -Inf))
            @test_throws "empty domains" ZigguratTools._check_arguments(256, (0, 0, 0))

            @test_throws "empty domains" ZigguratTools._check_arguments(256, [0, 0])
            @test_throws "empty domains" ZigguratTools._check_arguments(256, [Inf, Inf])
            @test_throws "empty domains" ZigguratTools._check_arguments(256, [-Inf, -Inf])
            @test_throws "empty domains" ZigguratTools._check_arguments(256, [0, 0, 0])
        end

        @testset "Infinite domains throw an infinite domain error" begin
            @test_throws "a domain" ZigguratTools._check_arguments(256, (-Inf, Inf))
            @test_throws "a domain" ZigguratTools._check_arguments(256, (Inf, -Inf))
            @test_throws "a domain" ZigguratTools._check_arguments(256, (Inf, 1, -Inf))
            @test_throws "a domain" ZigguratTools._check_arguments(256, (1, Inf, -Inf, 0))
            @test_throws "a domain" ZigguratTools._check_arguments(256, (-Inf, 1, Inf))

            @test_throws "a domain" ZigguratTools._check_arguments(256, [-Inf, Inf])
            @test_throws "a domain" ZigguratTools._check_arguments(256, [Inf, -Inf])
            @test_throws "a domain" ZigguratTools._check_arguments(256, [Inf, 1, -Inf])
            @test_throws "a domain" ZigguratTools._check_arguments(256, [1, Inf, -Inf, 0])
            @test_throws "a domain" ZigguratTools._check_arguments(256, [-Inf, 1, Inf])
        end

        @testset "Return nothing when everything checks out" begin
            @test ZigguratTools._check_arguments(256, (0, 1)) === nothing
            @test ZigguratTools._check_arguments(256, (0, 0.5, 1)) === nothing
            @test ZigguratTools._check_arguments(256, (1, 0, 0.5)) === nothing
            @test ZigguratTools._check_arguments(256, (1, 0)) === nothing
            @test ZigguratTools._check_arguments(256, (0, Inf)) === nothing
            @test ZigguratTools._check_arguments(256, (Inf, 0)) === nothing
            @test ZigguratTools._check_arguments(256, (0, 1, Inf)) === nothing
            @test ZigguratTools._check_arguments(256, (Inf, 1, 0)) === nothing
            @test ZigguratTools._check_arguments(256, (Inf, 0, 1)) === nothing
            @test ZigguratTools._check_arguments(256, (0, 0, 1)) === nothing
            @test ZigguratTools._check_arguments(256, (1, 0, 1)) === nothing
            @test ZigguratTools._check_arguments(256, (Inf, 0, Inf)) === nothing
            @test ZigguratTools._check_arguments(256, (0, Inf, 0)) === nothing
            @test ZigguratTools._check_arguments(256, (0, 1, Inf, Inf)) === nothing
            @test ZigguratTools._check_arguments(256, (Inf, 1, 0, 1)) === nothing
            @test ZigguratTools._check_arguments(256, (Inf, 0, 1, 1)) === nothing
        end
    end
end

@testset "_choose_tailarea_func" begin
    import ZigguratTools: _choose_tailarea_func
    inc_func(x) = x
    dec_func(x) = -x

    notnothing1 = x->1
    notnothing2 = x->2
    notnothing3 = x->3

    @test _choose_tailarea_func(inc_func, [0, 1], nothing, nothing, nothing) === nothing
    @test _choose_tailarea_func(dec_func, [0, 1], nothing, nothing, nothing) === nothing

    @test_throws "a ccdf is provided" _choose_tailarea_func(
        inc_func,
        [0, 1],
        nothing,
        nothing,
        notnothing1
    )
    @test _choose_tailarea_func(dec_func, [0, 1], nothing, nothing, notnothing1) ===
          notnothing1

    @test _choose_tailarea_func(inc_func, [0, 1], nothing, notnothing1, nothing) ===
          notnothing1
    @test_throws "a cdf is provided" _choose_tailarea_func(
        dec_func,
        [0, 1],
        nothing,
        notnothing1,
        nothing
    )

    @test _choose_tailarea_func(inc_func, [0, 1], nothing, notnothing1, notnothing2) ===
          notnothing1
    @test _choose_tailarea_func(dec_func, [0, 1], nothing, notnothing1, notnothing2) ===
          notnothing2

    @test _choose_tailarea_func(inc_func, [0, 1], notnothing1, nothing, nothing) ===
          notnothing1
    @test _choose_tailarea_func(dec_func, [0, 1], notnothing1, nothing, nothing) ===
          notnothing1

    @test _choose_tailarea_func(inc_func, [0, 1], notnothing1, nothing, notnothing2) ===
          notnothing1
    @test _choose_tailarea_func(dec_func, [0, 1], notnothing1, nothing, notnothing2) ===
          notnothing1

    @test _choose_tailarea_func(inc_func, [0, 1], notnothing1, notnothing2, nothing) ===
          notnothing1
    @test _choose_tailarea_func(dec_func, [0, 1], notnothing1, notnothing2, nothing) ===
          notnothing1

    @test _choose_tailarea_func(inc_func, [0, 1], notnothing1, notnothing2, notnothing3) ===
          notnothing1
    @test _choose_tailarea_func(dec_func, [0, 1], notnothing1, notnothing2, notnothing3) ===
          notnothing1
end
