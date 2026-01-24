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

            import Ziggurats: widths, layerratios, heights, highside
            @test widths(z_call) == widths(z_lamb) == widths(z_func)
            @test layerratios(z_call) == layerratios(z_lamb) == layerratios(z_func)
            @test heights(z_call) == heights(z_lamb) == heights(z_func)
            @test highside(z_call) == highside(z_lamb) == highside(z_func)
        end

        @testset "UnboundedZiggurat" begin
            z_call = UnboundedZiggurat(call, (0, Inf), 8)
            z_lamb = UnboundedZiggurat(lamb, (0, Inf), 8)
            z_func = UnboundedZiggurat(func, (0, Inf), 8)

            import Ziggurats: widths, layerratios, heights, highside
            @test widths(z_call) == widths(z_lamb) == widths(z_func)
            @test layerratios(z_call) == layerratios(z_lamb) == layerratios(z_func)
            @test heights(z_call) == heights(z_lamb) == heights(z_func)
            @test highside(z_call) == highside(z_lamb) == highside(z_func)
        end
    end

    @testset "ziggurat result depends on domain" begin
        f = x->exp(-x^2/2)

        @test ziggurat(f, (0, Inf), 2) isa UnboundedZiggurat{Float64}
        @test ziggurat(f, (0, Inf32), 2) isa UnboundedZiggurat{Float32}

        @test ziggurat(f, (0, 1), 2) isa BoundedZiggurat{Float64}
        @test ziggurat(f, (0, 1.0f0), 2) isa BoundedZiggurat{Float32}

        # Non-monotonic distributions should be detected.
        @test_throws "pdf is not monotonic" ziggurat(f, (-2, 1), 256)

        @test_broken try
            ziggurat(f, (-1, 1))
        catch e
            if startswith(e.msg, "pdf is not monotonic")
                true # test passes
            else
                rethrow()
            end
        else
            false # test fails
        end

        @test ziggurat(f, (-2, 0, 1), 2) isa CompositeZiggurat{Float64}
        @test_throws "pdf is not monotonic" ziggurat(f, (-2, 1, 1), 256)

        @test_broken try
            ziggurat(f, (-2, -1, 1), 2)
        catch e
            if startswith(e.msg, "pdf is not monotonic")
                true # test passes
            else
                rethrow()
            end
        else
            false # test fails
        end
    end

    @testset "ziggurat constructors fail fast when the domain is invalid" begin
        @testset "BoundedZiggurats fail when the domain contains Inf" begin
            invalid_domains = [(1, Inf), (-Inf, 1), (-Inf, Inf)]

            for d in invalid_domains
                @test_throws BoundedZiggurat(x -> 1, d, 4)
            end
        end

        @testset "UnboundedZiggurats fail when the domain does not contain exactly one Inf" begin
            invalid_domains = [(-Inf, Inf), (0, 1)]

            for d in invalid_domains
                @test_throws UnboundedZiggurat(x -> 1, d, 4)
            end
        end

        @testset "monotonic ziggurats fail when the domain is not a pair" begin
            invalid_domains = [(1, 2, 3), (1,), (1, 1.0)]

            constructors = [monotonic_ziggurat, BoundedZiggurat, UnboundedZiggurat]

            f = x->1

            @testset for d in invalid_domains, c in constructors
                @test_throws c(f, d, 4)
            end
        end

        @testset "domains without at least two distinct points are always invalid" begin
            @test_throws ziggurat(x->1, (), 4)
            @test_throws ziggurat(x->1, (1,), 4)
            @test_throws monotonic_ziggurat(x->1, (), 4)
            @test_throws monotonic_ziggurat(x->1, (1,), 4)
            @test_throws BellZiggurat(x->1, (), 4)
            @test_throws BellZiggurat(x->1, (1,), 4)
            @test_throws BoundedZiggurat(x->1, (), 4)
            @test_throws BoundedZiggurat(x->1, (1,), 4)
            @test_throws Unboundediggurat(x->1, (), 4)
            @test_throws Unboundediggurat(x->1, (1,), 4)
            @test_throws CompositeZiggurat(x->1, (), 4)
            @test_throws CompositeZiggurat(x->1, (1,), 4)
        end
    end
end
