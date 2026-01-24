issamevector(a::Vector{T}, b::Vector{T}) where {T} = a == b
issamevector(a::Vector, b::Vector) = false

@testset "Argument handling" begin
    @testset "regularize domains" begin
        import Ziggurats: regularize
        @testset "arrays, tuples, and iterators are accepted" begin
            domains = [[1, 2, 3], (1, 2, 3), 1:3, (i for i in 1:3)]

            @testset for d in domains
                @test issamevector(regularize(d).a, [1.0, 2.0, 3.0])
            end
        end

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
        using Ziggurats: _check_arguments
        valid_domains = regularize.([(0, Inf), (0, 1), (-Inf, 1)])

        inf_domains = regularize.([(-Inf, Inf)])

        wrong_length_domains = regularize.([(0, 0.5, 1), (0, 0.5, Inf), (-Inf, 0.5, 1), (-Inf, 0.5, Inf), ()])

        invalid_domains = [
            inf_domains;
            wrong_length_domains
        ]

        domains = [
            valid_domains;
            invalid_domains
        ]

        valid_Ns = [1, 2, 3, 4, 5, 2^6, 2^8, 2^12]
        invalid_Ns = [-1, 0]

        Ns = [
            valid_Ns;
            invalid_Ns
        ]

        @testset "Using a non-positive N throws an error" begin
            @testset for N in invalid_Ns, d in domains
                @test_throws DomainError _check_arguments(N, d)
            end
        end

        # Empty domian
        @testset "error when domain doesn't have exactly two distinct points" begin
            @testset for N in valid_Ns, d in wrong_length_domains
                @test_throws "the domian needs" _check_arguments(N, d)
            end
        end

        @testset "error when the domain is unbounded in both directions" begin
            @testset for N in valid_Ns, d in inf_domains
                @test_throws "a domain of (-Inf, Inf)" _check_arguments(N, d)
            end
        end

        @testset "Return nothing when everything checks out" begin
            for N in valid_Ns, d in valid_domains
                @test _check_arguments(N, d) === nothing
            end
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
