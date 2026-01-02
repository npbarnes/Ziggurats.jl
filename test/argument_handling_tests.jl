issamevector(a::Vector{T}, b::Vector{T}) where {T} = a == b
issamevector(a::Vector, b::Vector) = false

@testset "Argument handling" begin
    @testset "regularize domains" begin
        import Ziggurats: regularize
        @testset "regularize errors if it doesn't contain 2 or more distinct elements" begin
            @test_throws "empty domain" regularize(0)
            @test_throws "empty domain" regularize([0])
            @test_throws "empty domain" regularize((0,))
            @test_throws "empty domain" regularize([1, 1, 1])
            @test_throws "empty domain" regularize([1, 1.0f0, 1.0])
            @test_throws "empty domain" regularize([])
            @test_throws "empty domain" regularize(())
        end

        @testset "regularize promotes to highest float or Float64 if none" begin
            @test eltype(regularize((0, 1))) == Float64
            @test eltype(regularize((0, 1.0))) == Float64
            @test eltype(regularize((0.0, 1.0))) == Float64
            @test eltype(regularize((0.0f0, 1.0))) == Float64
            @test eltype(regularize((0.0, 1.0f0))) == Float64
            @test eltype(regularize((0, 1.0f0, 2.0))) == Float64

            @test eltype(regularize((0.0f0, 1))) == Float32
            @test eltype(regularize((0.0f0, 1.0f0))) == Float32
            @test eltype(regularize((0.0f0, Float16(1.0)))) == Float32
        end

        @testset "regularize sorts, promotes, and removes duplicates" begin
            @test issamevector(regularize((0, Inf)).a, [0.0, Inf])
            @test issamevector(regularize((Inf, 0)).a, [0.0, Inf])
            @test issamevector(regularize((0.0, Inf)).a, [0.0, Inf])
            @test issamevector(regularize((Inf, 0.0)).a, [0.0, Inf])
            @test issamevector(regularize((0.0f0, Inf)).a, [0.0, Inf])
            @test issamevector(regularize((Inf, 0.0f0)).a, [0.0, Inf])
            @test issamevector(regularize((0, 0.0, Inf)).a, [0.0, Inf])
            @test issamevector(regularize((0, 0.0, Inf32)).a, [0.0, Inf])
            @test issamevector(regularize((0, 0.0f0, Inf)).a, [0.0, Inf])
            @test issamevector(regularize((0, 0.0f0, Inf32)).a, [0.0f0, Inf32])
            @test issamevector(regularize((0.0f0, Inf32)).a, [0.0f0, Inf32])
            @test issamevector(regularize((Inf32, 0.0f0)).a, [0.0f0, Inf32])

            @test issamevector(regularize((0, 7)).a, [0.0, 7.0])
            @test issamevector(regularize((7, 0)).a, [0.0, 7.0])
            @test issamevector(regularize((0.0, 7)).a, [0.0, 7.0])
            @test issamevector(regularize((7.0, 0)).a, [0.0, 7.0])
            @test issamevector(regularize((0, 7.0)).a, [0.0, 7.0])
            @test issamevector(regularize((7, 0.0)).a, [0.0, 7.0])
            @test issamevector(regularize((0.0, 7.0)).a, [0.0, 7.0])
            @test issamevector(regularize((7.0, 0.0)).a, [0.0, 7.0])
            @test issamevector(regularize((0, 0, 7, 7)).a, [0.0, 7.0])
            @test issamevector(regularize((7, 7, 0, 0)).a, [0.0, 7.0])
            @test issamevector(regularize((0, 0.0, 7.0, 7)).a, [0.0, 7.0])
            @test issamevector(regularize((7, 7.0, 0.0, 0)).a, [0.0, 7.0])

            @test issamevector(regularize((0.0f0, 7)).a, [0.0f0, 7.0f0])
            @test issamevector(regularize((7.0f0, 0)).a, [0.0f0, 7.0f0])
            @test issamevector(regularize((0, 7.0f0)).a, [0.0f0, 7.0f0])
            @test issamevector(regularize((7, 0.0f0)).a, [0.0f0, 7.0f0])
            @test issamevector(regularize((0.0f0, 7.0f0)).a, [0.0f0, 7.0f0])
            @test issamevector(regularize((7.0f0, 0.0f0)).a, [0.0f0, 7.0f0])
            @test issamevector(regularize((0, 0.0f0, 7.0f0, 7)).a, [0.0f0, 7.0f0])
            @test issamevector(regularize((7, 7.0f0, 0.0f0, 0)).a, [0.0f0, 7.0f0])
        end

        @testset "regularize is idempotent" begin
            r = regularize((1, 2, 3))
            @test r === regularize(r)

            r = regularize((1.0, 2, 3))
            @test r === regularize(r)

            r = regularize((1.0f0, 2, 3))
            @test r === regularize(r)

            r = regularize([1, 2, 3])
            @test r === regularize(r)

            r = regularize([1.0, 2, 3])
            @test r === regularize(r)

            r = regularize([1.0f0, 2, 3])
            @test r === regularize(r)
        end

        @testset "regularize copies an array" begin
            a = [1.0, 2.0, 3.0]
            ra = regularize(a).a

            @test issamevector(a, ra)
            @test a !== ra
        end
    end

    @testset "_check_arguments" begin
        @testset "Using a non-positive N throws an error" begin
            @testset "When domain is valid" begin
                @test_throws "N must" Ziggurats._check_arguments(-1, (0, Inf))
                @test_throws "N must" Ziggurats._check_arguments(-1, (0, 1))
                @test_throws "N must" Ziggurats._check_arguments(-1, (-Inf, 1))
                @test_throws "N must" Ziggurats._check_arguments(-1, (-Inf, Inf))

                @test_throws "N must" Ziggurats._check_arguments(-1, [0, 1])
                @test_throws "N must" Ziggurats._check_arguments(-1, [0, Inf])
                @test_throws "N must" Ziggurats._check_arguments(-1, [-Inf, 1])
                @test_throws "N must" Ziggurats._check_arguments(-1, [-Inf, Inf])

                @test_throws "N must" Ziggurats._check_arguments(0, (0, 1))
                @test_throws "N must" Ziggurats._check_arguments(0, (0, Inf))
                @test_throws "N must" Ziggurats._check_arguments(0, (-Inf, 1))
                @test_throws "N must" Ziggurats._check_arguments(0, (-Inf, Inf))

                @test_throws "N must" Ziggurats._check_arguments(0, [0, 1])
                @test_throws "N must" Ziggurats._check_arguments(0, [0, Inf])
                @test_throws "N must" Ziggurats._check_arguments(0, [-Inf, 1])
                @test_throws "N must" Ziggurats._check_arguments(0, [-Inf, Inf])

                @test_throws "N must" Ziggurats._check_arguments(-1, (0, 0.5, 1))
                @test_throws "N must" Ziggurats._check_arguments(-1, (0, 0.5, Inf))
                @test_throws "N must" Ziggurats._check_arguments(-1, (-Inf, 0.5, 1))
                @test_throws "N must" Ziggurats._check_arguments(-1, (-Inf, 0.5, Inf))

                @test_throws "N must" Ziggurats._check_arguments(-1, [0, 0.5, 1])
                @test_throws "N must" Ziggurats._check_arguments(-1, [0, 0.5, Inf])
                @test_throws "N must" Ziggurats._check_arguments(-1, [-Inf, 0.5, 1])
                @test_throws "N must" Ziggurats._check_arguments(-1, [-Inf, 0.5, Inf])

                @test_throws "N must" Ziggurats._check_arguments(0, (0, 0.5, 1))
                @test_throws "N must" Ziggurats._check_arguments(0, (0, 0.5, Inf))
                @test_throws "N must" Ziggurats._check_arguments(0, (-Inf, 0.5, 1))
                @test_throws "N must" Ziggurats._check_arguments(0, (-Inf, 0.5, Inf))

                @test_throws "N must" Ziggurats._check_arguments(0, [0, 0.5, 1])
                @test_throws "N must" Ziggurats._check_arguments(0, [0, 0.5, Inf])
                @test_throws "N must" Ziggurats._check_arguments(0, [-Inf, 0.5, 1])
                @test_throws "N must" Ziggurats._check_arguments(0, [-Inf, 0.5, Inf])
            end

            @testset "When the domain is invalid" begin
                @test_throws "N must be" Ziggurats._check_arguments(-1, ())
                @test_throws "N must be" Ziggurats._check_arguments(-1, (0,))
                @test_throws "N must be" Ziggurats._check_arguments(-1, (Inf,))
                @test_throws "N must be" Ziggurats._check_arguments(-1, (0, 0))
                @test_throws "N must be" Ziggurats._check_arguments(-1, (Inf, Inf))
                @test_throws "N must be" Ziggurats._check_arguments(-1, (-Inf, -Inf))
                @test_throws "N must be" Ziggurats._check_arguments(-1, (-Inf, Inf))

                @test_throws "N must be" Ziggurats._check_arguments(0, (0, 0))
                @test_throws "N must be" Ziggurats._check_arguments(0, (Inf, Inf))
                @test_throws "N must be" Ziggurats._check_arguments(0, (-Inf, -Inf))
                @test_throws "N must be" Ziggurats._check_arguments(0, (-Inf, Inf))

                @test_throws "N must be" Ziggurats._check_arguments(-1, [])
                @test_throws "N must be" Ziggurats._check_arguments(-1, [0])
                @test_throws "N must be" Ziggurats._check_arguments(-1, [Inf])
                @test_throws "N must be" Ziggurats._check_arguments(-1, [0, 0])
                @test_throws "N must be" Ziggurats._check_arguments(-1, [Inf, Inf])
                @test_throws "N must be" Ziggurats._check_arguments(-1, [-Inf, -Inf])
                @test_throws "N must be" Ziggurats._check_arguments(-1, [-Inf, Inf])

                @test_throws "N must be" Ziggurats._check_arguments(0, [0, 0])
                @test_throws "N must be" Ziggurats._check_arguments(0, [Inf, Inf])
                @test_throws "N must be" Ziggurats._check_arguments(0, [-Inf, -Inf])
                @test_throws "N must be" Ziggurats._check_arguments(0, [-Inf, Inf])
            end
        end

        # Empty domian
        @testset "Empty domains throw an empty domain error" begin
            @test_throws "empty domains" Ziggurats._check_arguments(256, ())
            @test_throws "empty domains" Ziggurats._check_arguments(256, (0,))
            @test_throws "empty domains" Ziggurats._check_arguments(256, (Inf,))
            @test_throws "empty domains" Ziggurats._check_arguments(256, (0, 0))
            @test_throws "empty domains" Ziggurats._check_arguments(256, (Inf, Inf))
            @test_throws "empty domains" Ziggurats._check_arguments(256, (-Inf, -Inf))
            @test_throws "empty domains" Ziggurats._check_arguments(256, (0, 0, 0))

            @test_throws "empty domains" Ziggurats._check_arguments(256, [])
            @test_throws "empty domains" Ziggurats._check_arguments(256, [0])
            @test_throws "empty domains" Ziggurats._check_arguments(256, [Inf])
            @test_throws "empty domains" Ziggurats._check_arguments(256, [0, 0])
            @test_throws "empty domains" Ziggurats._check_arguments(256, [Inf, Inf])
            @test_throws "empty domains" Ziggurats._check_arguments(256, [-Inf, -Inf])
            @test_throws "empty domains" Ziggurats._check_arguments(256, [0, 0, 0])
        end

        @testset "Infinite domains throw an infinite domain error" begin
            @test_throws "a domain" Ziggurats._check_arguments(256, (-Inf, Inf))
            @test_throws "a domain" Ziggurats._check_arguments(256, (Inf, -Inf))
            @test_throws "a domain" Ziggurats._check_arguments(256, (Inf, 1, -Inf))
            @test_throws "a domain" Ziggurats._check_arguments(256, (1, Inf, -Inf, 0))
            @test_throws "a domain" Ziggurats._check_arguments(256, (-Inf, 1, Inf))

            @test_throws "a domain" Ziggurats._check_arguments(256, [-Inf, Inf])
            @test_throws "a domain" Ziggurats._check_arguments(256, [Inf, -Inf])
            @test_throws "a domain" Ziggurats._check_arguments(256, [Inf, 1, -Inf])
            @test_throws "a domain" Ziggurats._check_arguments(256, [1, Inf, -Inf, 0])
            @test_throws "a domain" Ziggurats._check_arguments(256, [-Inf, 1, Inf])
        end

        @testset "Return nothing when everything checks out" begin
            @test Ziggurats._check_arguments(256, (0, 1)) === nothing
            @test Ziggurats._check_arguments(256, (0, 0.5, 1)) === nothing
            @test Ziggurats._check_arguments(256, (1, 0, 0.5)) === nothing
            @test Ziggurats._check_arguments(256, (1, 0)) === nothing
            @test Ziggurats._check_arguments(256, (0, Inf)) === nothing
            @test Ziggurats._check_arguments(256, (Inf, 0)) === nothing
            @test Ziggurats._check_arguments(256, (0, 1, Inf)) === nothing
            @test Ziggurats._check_arguments(256, (Inf, 1, 0)) === nothing
            @test Ziggurats._check_arguments(256, (Inf, 0, 1)) === nothing
            @test Ziggurats._check_arguments(256, (0, 0, 1)) === nothing
            @test Ziggurats._check_arguments(256, (1, 0, 1)) === nothing
            @test Ziggurats._check_arguments(256, (Inf, 0, Inf)) === nothing
            @test Ziggurats._check_arguments(256, (0, Inf, 0)) === nothing
            @test Ziggurats._check_arguments(256, (0, 1, Inf, Inf)) === nothing
            @test Ziggurats._check_arguments(256, (Inf, 1, 0, 1)) === nothing
            @test Ziggurats._check_arguments(256, (Inf, 0, 1, 1)) === nothing
        end
    end

    @testset "_choose_tailarea_func" begin
        import Ziggurats: _choose_tailarea_func
        inc_func(x) = x
        dec_func(x) = -x

        notnothing1 = x->1
        notnothing2 = x->2
        notnothing3 = x->3

        @test _choose_tailarea_func(inc_func, [0, 1], nothing, nothing, nothing) === nothing
        @test _choose_tailarea_func(dec_func, [0, 1], nothing, nothing, nothing) === nothing

        @test_throws "a ccdf is provided" _choose_tailarea_func(inc_func, [0, 1], nothing, nothing, notnothing1)
        @test _choose_tailarea_func(dec_func, [0, 1], nothing, nothing, notnothing1) === notnothing1

        @test _choose_tailarea_func(inc_func, [0, 1], nothing, notnothing1, nothing) === notnothing1
        @test_throws "a cdf is provided" _choose_tailarea_func(dec_func, [0, 1], nothing, notnothing1, nothing)

        @test _choose_tailarea_func(inc_func, [0, 1], nothing, notnothing1, notnothing2) === notnothing1
        @test _choose_tailarea_func(dec_func, [0, 1], nothing, notnothing1, notnothing2) === notnothing2

        @test _choose_tailarea_func(inc_func, [0, 1], notnothing1, nothing, nothing) === notnothing1
        @test _choose_tailarea_func(dec_func, [0, 1], notnothing1, nothing, nothing) === notnothing1

        @test _choose_tailarea_func(inc_func, [0, 1], notnothing1, nothing, notnothing2) === notnothing1
        @test _choose_tailarea_func(dec_func, [0, 1], notnothing1, nothing, notnothing2) === notnothing1

        @test _choose_tailarea_func(inc_func, [0, 1], notnothing1, notnothing2, nothing) === notnothing1
        @test _choose_tailarea_func(dec_func, [0, 1], notnothing1, notnothing2, nothing) === notnothing1

        @test _choose_tailarea_func(inc_func, [0, 1], notnothing1, notnothing2, notnothing3) === notnothing1
        @test _choose_tailarea_func(dec_func, [0, 1], notnothing1, notnothing2, notnothing3) === notnothing1
    end

    @testset "_identify_mode" begin
        import Ziggurats: _identify_mode
        f = x->x
        @test _identify_mode(f, (0, 1)) == (1, 0)

        f = x->5-x
        @test _identify_mode(f, (0, 1)) == (0, 1)

        @testset "the finite side of an unbounded domain is the mode" begin
            f = x->error("This function shouldn't be evaluated")
            @test _identify_mode(f, (0, Inf)) == (0, Inf)
            @test _identify_mode(f, (-Inf, 0)) == (0, -Inf)
        end

        @testset "Check for invalid values on the boundaries" begin
            f = x -> x==1 ? NaN : x
            @test_throws "pdf(x) is NaN on the boundary." _identify_mode(f, (0, 1))
            @test_throws "pdf(x) is NaN on the boundary." _identify_mode(f, (1, 2))

            f = x -> x==1 ? Inf : x
            @test_throws "pdf(x) is infinite on the boundary." _identify_mode(f, (0, 1))
            @test_throws "pdf(x) is infinite on the boundary." _identify_mode(f, (1, 2))
        end
    end
end
