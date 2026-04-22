@testset "Interface Tests" begin
    @testset "Different types of callables are accepted" begin
        struct Callable end
        (c::Callable)(x) = exp(-x^2)

        call = Callable()
        lamb = x -> exp(-x^2)
        func(x) = exp(-x^2)

        @testset "BoundedZiggurat" begin
            z_call = BoundedZiggurat{Float64,Float64,8}(call, (0, 1))
            z_lamb = BoundedZiggurat{Float64,Float64,8}(lamb, (0, 1))
            z_func = BoundedZiggurat{Float64,Float64,8}(func, (0, 1))

            import Ziggurats: widths, layerratios, heights, highside
            @test widths(z_call) == widths(z_lamb) == widths(z_func)
            @test layerratios(z_call) == layerratios(z_lamb) == layerratios(z_func)
            @test heights(z_call) == heights(z_lamb) == heights(z_func)
            @test highside(z_call) == highside(z_lamb) == highside(z_func)
        end

        @testset "UnboundedZiggurat" begin
            z_call = UnboundedZiggurat{Float64,Float64,8}(call, (0, Inf))
            z_lamb = UnboundedZiggurat{Float64,Float64,8}(lamb, (0, Inf))
            z_func = UnboundedZiggurat{Float64,Float64,8}(func, (0, Inf))

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

        @test ziggurat(f, (-2, 0, 1), 2) isa Ziggurats.CompositeZiggurat{Float64}
        @test_throws "pdf is not monotonic" ziggurat(f, (-2, 1), 256)

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
            @test_throws "expected a bounded domain" BoundedZiggurat{Float64,Float64,4}(x -> 1, (1, Inf))
            @test_throws "expected a bounded domain" BoundedZiggurat{Float64,Float64,4}(x -> 1, (-Inf, 1))
            @test_throws "expected a bounded domain" BoundedZiggurat{Float64,Float64,4}(x -> 1, (1, Inf))
            @test_throws "a domain of (-Inf, Inf) is impossible" BoundedZiggurat{Float64,Float64,4}(x -> 1, (-Inf, Inf))
        end

        @testset "UnboundedZiggurats fail when the domain does not contain exactly one Inf" begin
            @test_throws "expected an unbounded domain" UnboundedZiggurat{Float64,Float64,4}(x -> 1, (0, 1))
            @test_throws "a domain of (-Inf, Inf) is impossible" UnboundedZiggurat{Float64,Float64,4}(
                x -> 1,
                (-Inf, Inf)
            )
        end

        @testset "monotonic ziggurats fail when the domain is not a pair" begin
            f = x->1

            @test_throws ErrorException monotonic_ziggurat(f, (1, 2, 3), 4)
            @test_throws ErrorException BoundedZiggurat{Float64,Float64,4}(f, (1, 2, 3))
            @test_throws ErrorException UnboundedZiggurat{Float64,Float64,4}(f, (1, 2, Inf))
        end

        @testset "domains without at least two distinct points are always invalid" begin
            @test_throws ErrorException ziggurat(x->1, (), 4)
            @test_throws ErrorException ziggurat(x->1, (1,), 4)
            @test_throws ErrorException monotonic_ziggurat(x->1, (), 4)
            @test_throws ErrorException monotonic_ziggurat(x->1, (1,), 4)
            @test_throws ErrorException BoundedBellZiggurat{Float64,Float64,4}(x->1, ())
            @test_throws ErrorException BoundedBellZiggurat{Float64,Float64,4}(x->1, (1,))
            @test_throws ErrorException UnboundedBellZiggurat{Float64,Float64,4}(x->1, ())
            @test_throws ErrorException UnboundedBellZiggurat{Float64,Float64,4}(x->1, (1,))
            @test_throws ErrorException BoundedZiggurat{Float64,Float64,4}(x->1, ())
            @test_throws ErrorException BoundedZiggurat{Float64,Float64,4}(x->1, (1,))
            @test_throws ErrorException UnboundedZiggurat{Float64,Float64,4}(x->1, ())
            @test_throws ErrorException UnboundedZiggurat{Float64,Float64,4}(x->1, (1,))
        end
    end
end
