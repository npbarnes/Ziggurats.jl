@testset "Type-stable interface" begin
    @testset "Capital letter constructors are always type stable" begin
        using Ziggurats: default_numlayers
        # Use functions that are not type stable or inferrable and domains whose entries are
        # not type inferrable or stable as a worst case scenario for type-stability
        unstable_f = x -> begin
            r = exp(-abs(x))
            r = reinterpret(Unsigned, x) % Bool ? BigFloat(r) : Float64(r)
            Ref{Any}(r)[]
        end
        unstable_ipdf_left = y -> begin
            r = log(y)
            r = reinterpret(Unsigned, y) % Bool ? BigFloat(r) : Float64(r)
            Ref{Any}(r)[]
        end
        unstable_ipdf_right = y -> begin
            r = -log(y)
            r = reinterpret(Unsigned, y) % Bool ? BigFloat(r) : Float64(r)
            Ref{Any}(r)[]
        end
        unstable_cdf = x -> begin
            if x < zero(x)
                r = exp(x)/2
            else
                r = one(x) - exp(-x)/2
            end
            r = reinterpret(Unsigned, x) % Bool ? BigFloat(r) : Float64(r)
            Ref{Any}(r)[]
        end
        unstable_ccdf = x -> begin
            if x < zero(x)
                r = one(x) - exp(x)/2
            else
                r = exp(-x)/2
            end
            r = reinterpret(Unsigned, x) % Bool ? BigFloat(r) : Float64(r)
            Ref{Any}(r)[]
        end
        unstable_fallback_left = (rng, a) -> begin
            u = rand(rng, typeof(a))
            r = a + log1p(-u)
            r = reinterpret(Unsigned, u) % Bool ? BigFloat(r) : Float64(r)
            Ref{Any}(r)[]
        end
        unstable_fallback_right = (rng, a) -> begin
            u = rand(rng, typeof(a))
            r = a - log1p(-u)
            r = reinterpret(Unsigned, u) % Bool ? BigFloat(r) : Float64(r)
            Ref{Any}(r)[]
        end
        @testset "BoundedZiggurat" begin
            @testset for X in TestTypes, Y in TestTypes
                @testset for ipdf in (nothing, unstable_ipdf_right)
                    # Default num layers
                    @test @inferred(BoundedZiggurat{X,Y}(unstable_f, Any[0, 1.0]; ipdf)) isa
                          BoundedZiggurat{X,Y,default_numlayers(X, Y)}

                    @test @inferred(BoundedZiggurat{X,Y,nothing}(unstable_f, Any[0, 1.0]; ipdf)) isa
                          BoundedZiggurat{X,Y,default_numlayers(X, Y)}

                    # Power of two
                    @test @inferred(BoundedZiggurat{X,Y,4}(unstable_f, Any[0, 1.0]; ipdf)) isa BoundedZiggurat{X,Y,4}

                    # Non-power of two
                    @test @inferred(BoundedZiggurat{X,Y,5}(unstable_f, Any[0, 1.0]; ipdf)) isa BoundedZiggurat{X,Y,5}
                end
            end
        end
        @testset "UnboundedZiggurat" begin
            @testset for X in TestTypes, Y in TestTypes
                @testset for (ipdf, tailarea, fallback) in (
                    (nothing, nothing, nothing),
                    (unstable_ipdf_right, nothing, nothing),
                    (nothing, unstable_ccdf, nothing),
                    (nothing, nothing, unstable_fallback_right)
                )
                    # Default num layers
                    @test @inferred(UnboundedZiggurat{X,Y}(unstable_f, Any[0, Inf]; ipdf, tailarea, fallback)) isa
                          UnboundedZiggurat{X,Y,default_numlayers(X, Y)}

                    @test @inferred(UnboundedZiggurat{X,Y,nothing}(
                        unstable_f,
                        Any[0, Inf];
                        ipdf,
                        tailarea,
                        fallback
                    )) isa UnboundedZiggurat{X,Y,default_numlayers(X, Y)}

                    # Power of two
                    @test @inferred(UnboundedZiggurat{X,Y,4}(unstable_f, Any[0, Inf]; ipdf, tailarea, fallback)) isa
                          UnboundedZiggurat{X,Y,4}

                    # Non-power of two
                    @test @inferred(UnboundedZiggurat{X,Y,5}(unstable_f, Any[0, Inf]; ipdf, tailarea, fallback)) isa
                          UnboundedZiggurat{X,Y,5}
                end
            end
        end
        @testset "BoundedBellZiggurat" begin
            @testset for X in TestTypes, Y in TestTypes
                @testset for ipdf in (nothing, unstable_ipdf_right)
                    # Default num layers
                    @test @inferred(BoundedBellZiggurat{X,Y}(unstable_f, Any[0, 1.0]; ipdf)) isa
                          BoundedBellZiggurat{X,Y,default_numlayers(X, Y)}

                    @test @inferred(BoundedBellZiggurat{X,Y,nothing}(unstable_f, Any[0, 1.0]; ipdf)) isa
                          BoundedBellZiggurat{X,Y,default_numlayers(X, Y)}

                    # Power of two
                    @test @inferred(BoundedBellZiggurat{X,Y,4}(unstable_f, Any[0, 1.0]; ipdf)) isa
                          BoundedBellZiggurat{X,Y,4}

                    # Non-power of two
                    @test @inferred(BoundedBellZiggurat{X,Y,5}(unstable_f, Any[0, 1.0]; ipdf)) isa
                          BoundedBellZiggurat{X,Y,5}
                end
            end
        end
        @testset "UnboundedBellZiggurat" begin
            @testset for X in TestTypes, Y in TestTypes
                @testset for (ipdf, tailarea, fallback) in (
                    (nothing, nothing, nothing),
                    (unstable_ipdf_right, nothing, nothing),
                    (nothing, unstable_ccdf, nothing),
                    (nothing, nothing, unstable_fallback_right)
                )
                    # Default num layers
                    @test @inferred(UnboundedBellZiggurat{X,Y}(unstable_f, Any[0, Inf]; ipdf, tailarea, fallback)) isa
                          UnboundedBellZiggurat{X,Y,default_numlayers(X, Y)}

                    @test @inferred(UnboundedBellZiggurat{X,Y,nothing}(
                        unstable_f,
                        Any[0, Inf];
                        ipdf,
                        tailarea,
                        fallback
                    )) isa UnboundedBellZiggurat{X,Y,default_numlayers(X, Y)}

                    # Power of two
                    @test @inferred(UnboundedBellZiggurat{X,Y,4}(unstable_f, Any[0, Inf]; ipdf, tailarea, fallback)) isa
                          UnboundedBellZiggurat{X,Y,4}

                    # Non-power of two
                    @test @inferred(UnboundedBellZiggurat{X,Y,5}(unstable_f, Any[0, Inf]; ipdf, tailarea, fallback)) isa
                          UnboundedBellZiggurat{X,Y,5}
                end
            end
        end
        @testset "BoundedCompositeZiggurat" begin
            @testset "Errors" begin
                @test_throws MethodError BoundedCompositeZiggurat{Float32,Float16}(unstable_f, [-1, 0, 1])
                @test_throws MethodError BoundedCompositeZiggurat{Float32,Float16}(unstable_f, Any[-1, 0, 1])
                @test_throws MethodError BoundedCompositeZiggurat{Float32,Float16,nothing}(unstable_f, Any[-1, 0, 1])
                @test_throws MethodError BoundedCompositeZiggurat{Float32,Float16,3}(unstable_f, Any[-1, 0, 1])
                @test_throws MethodError BoundedCompositeZiggurat{Float32,Float16,4}(unstable_f, Any[-1, 0, 1])
                @test_throws MethodError BoundedCompositeZiggurat{Float32,Float16,(nothing, nothing)}(
                    unstable_f,
                    Any[-1, 0, 1]
                )

                @test_throws ErrorException BoundedCompositeZiggurat{Float32,Float16,Tuple{nothing,nothing},3}(
                    unstable_f,
                    Any[-1, 0, 1]
                )
                @test_throws ErrorException BoundedCompositeZiggurat{Float32,Float16,Tuple{nothing,nothing},1}(
                    unstable_f,
                    Any[-1, 0, 1]
                )
                @test_throws ErrorException BoundedCompositeZiggurat{Float32,Float16,NTuple{4,nothing},5}(
                    unstable_f,
                    Any[-1, 0, 0.5, 1, 1.5]
                )
                @test_throws ErrorException BoundedCompositeZiggurat{Float32,Float16,NTuple{4,nothing},3}(
                    unstable_f,
                    Any[-1, 0, 0.5, 1, 1.5]
                )
            end
            @testset for X in (Float64, Float32), Y in (Float64, Float32)
                @testset for (ipdf_left, ipdf_right) in ((nothing, nothing), (unstable_ipdf_left, unstable_ipdf_right)),
                    cdf in (nothing, unstable_cdf)
                    # Two parameter forms
                    @test @inferred(BoundedCompositeZiggurat{X,Y}(
                        unstable_f,
                        (-1, 0, 1);
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,NTuple{2,default_numlayers(X, Y)}}

                    @test @inferred(BoundedCompositeZiggurat{X,Y}(
                        unstable_f,
                        SVector{3,Any}(-1, 0, 1);
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,NTuple{2,default_numlayers(X, Y)}}

                    @test @inferred(BoundedCompositeZiggurat{X,Y}(
                        unstable_f,
                        (-1, 0, 0.5, 1);
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,NTuple{3,default_numlayers(X, Y)}}

                    @test @inferred(BoundedCompositeZiggurat{X,Y}(
                        unstable_f,
                        SVector{4,Any}(-1, 0, 0.5, 1);
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,NTuple{3,default_numlayers(X, Y)}}

                    # Three parameter forms
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,nothing}}(
                        unstable_f,
                        Any[-1, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{3,nothing}}(
                        unstable_f,
                        Any[-1, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,Tuple{3,default_numlayers(X, Y)}}
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,3}}(
                        unstable_f,
                        Any[-1, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),3}}

                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,nothing}}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{
                        X,
                        Y,
                        Tuple{
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y)
                        }
                    }
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{3,nothing,nothing,nothing}}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{
                        X,
                        Y,
                        Tuple{3,default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,3,nothing,nothing}}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),3,default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,3}}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y),3}
                    }

                    # Four parameter forms
                    @test @inferred(BoundedCompositeZiggurat{X,Y,nothing,2}(
                        unstable_f,
                        Any[-1, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}

                    @test @inferred(BoundedCompositeZiggurat{X,Y,nothing,4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,NTuple{4,default_numlayers(X, Y)}}

                    @test @inferred(BoundedCompositeZiggurat{X,Y,3,2}(
                        unstable_f,
                        Any[-1, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,Tuple{3,3}}

                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,nothing},2}(
                        unstable_f,
                        Any[-1, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{3,nothing},2}(
                        unstable_f,
                        Any[-1, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{X,Y,Tuple{3,default_numlayers(X, Y)}}

                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,nothing},4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{
                        X,
                        Y,
                        Tuple{
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y)
                        }
                    }
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{3,nothing,nothing,nothing},4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{
                        X,
                        Y,
                        Tuple{3,default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,3,nothing,nothing},4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),3,default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(BoundedCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,3},4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf
                    )) isa BoundedCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y),3}
                    }
                end
            end
        end
        @testset "LeftTailCompositeZiggurat" begin
            @testset "Errors" begin
                @test_throws MethodError LeftTailCompositeZiggurat{Float32,Float16}(unstable_f, [-Inf, 0, 1])
                @test_throws MethodError LeftTailCompositeZiggurat{Float32,Float16}(unstable_f, Any[-Inf, 0, 1])
                @test_throws MethodError LeftTailCompositeZiggurat{Float32,Float16,nothing}(unstable_f, Any[-Inf, 0, 1])
                @test_throws MethodError LeftTailCompositeZiggurat{Float32,Float16,3}(unstable_f, Any[-Inf, 0, 1])
                @test_throws MethodError LeftTailCompositeZiggurat{Float32,Float16,4}(unstable_f, Any[-Inf, 0, 1])
                @test_throws MethodError LeftTailCompositeZiggurat{Float32,Float16,(nothing, nothing)}(
                    unstable_f,
                    Any[-Inf, 0, 1]
                )

                @test_throws ErrorException LeftTailCompositeZiggurat{Float32,Float16,Tuple{nothing,nothing},3}(
                    unstable_f,
                    Any[-Inf, 0, 1]
                )
                @test_throws ErrorException LeftTailCompositeZiggurat{Float32,Float16,Tuple{nothing,nothing},1}(
                    unstable_f,
                    Any[-Inf, 0, 1]
                )
                @test_throws ErrorException LeftTailCompositeZiggurat{Float32,Float16,NTuple{4,nothing},5}(
                    unstable_f,
                    Any[-Inf, 0, 0.5, 1, 1.5]
                )
                @test_throws ErrorException LeftTailCompositeZiggurat{Float32,Float16,NTuple{4,nothing},3}(
                    unstable_f,
                    Any[-Inf, 0, 0.5, 1, 1.5]
                )
            end
            @testset for X in (Float64, Float32), Y in (Float64, Float32)
                @testset for (ipdf_left, ipdf_right) in ((nothing, nothing), (unstable_ipdf_left, unstable_ipdf_right)),
                    cdf in (nothing, unstable_cdf),
                    left_fallback in (nothing, unstable_fallback_left)
                    # Two parameter forms
                    @test @inferred(LeftTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        (-Inf, 0, 1);
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,NTuple{2,default_numlayers(X, Y)}}

                    @test @inferred(LeftTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        SVector{3,Any}(-Inf, 0, 1);
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,NTuple{2,default_numlayers(X, Y)}}

                    @test @inferred(LeftTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        (-Inf, 0, 0.5, 1);
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,NTuple{3,default_numlayers(X, Y)}}

                    @test @inferred(LeftTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        SVector{4,Any}(-Inf, 0, 0.5, 1);
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,NTuple{3,default_numlayers(X, Y)}}

                    # Three parameter forms
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{3,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,Tuple{3,default_numlayers(X, Y)}}
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,3}}(
                        unstable_f,
                        Any[-Inf, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),3}}

                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y)
                        }
                    }
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{3,nothing,nothing,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{3,default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,3,nothing,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),3,default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,3}}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y),3}
                    }

                    # Four parameter forms
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,nothing,2}(
                        unstable_f,
                        Any[-Inf, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}

                    @test @inferred(LeftTailCompositeZiggurat{X,Y,nothing,4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,NTuple{4,default_numlayers(X, Y)}}

                    @test @inferred(LeftTailCompositeZiggurat{X,Y,3,2}(
                        unstable_f,
                        Any[-Inf, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,Tuple{3,3}}

                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,nothing},2}(
                        unstable_f,
                        Any[-Inf, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{3,nothing},2}(
                        unstable_f,
                        Any[-Inf, 0, 1];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{X,Y,Tuple{3,default_numlayers(X, Y)}}

                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,nothing},4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y)
                        }
                    }
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{3,nothing,nothing,nothing},4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{3,default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,3,nothing,nothing},4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),3,default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(LeftTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,3},4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, 1.5];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        left_fallback
                    )) isa LeftTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y),3}
                    }
                end
            end
        end
        @testset "RightTailCompositeZiggurat" begin
            @testset "Errors" begin
                @test_throws MethodError RightTailCompositeZiggurat{Float32,Float16}(unstable_f, [-1, 0, Inf])
                @test_throws MethodError RightTailCompositeZiggurat{Float32,Float16}(unstable_f, Any[-1, 0, Inf])
                @test_throws MethodError RightTailCompositeZiggurat{Float32,Float16,nothing}(
                    unstable_f,
                    Any[-1, 0, Inf]
                )
                @test_throws MethodError RightTailCompositeZiggurat{Float32,Float16,3}(unstable_f, Any[-1, 0, Inf])
                @test_throws MethodError RightTailCompositeZiggurat{Float32,Float16,4}(unstable_f, Any[-1, 0, Inf])
                @test_throws MethodError RightTailCompositeZiggurat{Float32,Float16,(nothing, nothing)}(
                    unstable_f,
                    Any[-1, 0, Inf]
                )

                @test_throws ErrorException RightTailCompositeZiggurat{Float32,Float16,Tuple{nothing,nothing},3}(
                    unstable_f,
                    Any[-1, 0, Inf]
                )
                @test_throws ErrorException RightTailCompositeZiggurat{Float32,Float16,Tuple{nothing,nothing},1}(
                    unstable_f,
                    Any[-1, 0, Inf]
                )
                @test_throws ErrorException RightTailCompositeZiggurat{Float32,Float16,NTuple{4,nothing},5}(
                    unstable_f,
                    Any[-1, 0, 0.5, 1, Inf]
                )
                @test_throws ErrorException RightTailCompositeZiggurat{Float32,Float16,NTuple{4,nothing},3}(
                    unstable_f,
                    Any[-1, 0, 0.5, 1, Inf]
                )
            end
            @testset for X in (Float64, Float32), Y in (Float64, Float32)
                @testset for (ipdf_left, ipdf_right) in ((nothing, nothing), (unstable_ipdf_left, unstable_ipdf_right)),
                    ccdf in (nothing, unstable_ccdf),
                    right_fallback in (nothing, unstable_fallback_right)
                    # Two parameter forms
                    @test @inferred(RightTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        (-1, 0, Inf);
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,NTuple{2,default_numlayers(X, Y)}}

                    @test @inferred(RightTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        SVector{3,Any}(-1, 0, Inf);
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,NTuple{2,default_numlayers(X, Y)}}

                    @test @inferred(RightTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        (-1, 0, 0.5, Inf);
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,NTuple{3,default_numlayers(X, Y)}}

                    @test @inferred(RightTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        SVector{4,Any}(-1, 0, 0.5, Inf);
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,NTuple{3,default_numlayers(X, Y)}}

                    # Three parameter forms
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,nothing}}(
                        unstable_f,
                        Any[-1, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{3,nothing}}(
                        unstable_f,
                        Any[-1, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,Tuple{3,default_numlayers(X, Y)}}
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,3}}(
                        unstable_f,
                        Any[-1, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),3}}

                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,nothing}}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y)
                        }
                    }
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{3,nothing,nothing,nothing}}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{3,default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,3,nothing,nothing}}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),3,default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,3}}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y),3}
                    }

                    # Four parameter forms
                    @test @inferred(RightTailCompositeZiggurat{X,Y,nothing,2}(
                        unstable_f,
                        Any[-1, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}

                    @test @inferred(RightTailCompositeZiggurat{X,Y,nothing,4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,NTuple{4,default_numlayers(X, Y)}}

                    @test @inferred(RightTailCompositeZiggurat{X,Y,3,2}(
                        unstable_f,
                        Any[-1, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,Tuple{3,3}}

                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,nothing},2}(
                        unstable_f,
                        Any[-1, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{3,nothing},2}(
                        unstable_f,
                        Any[-1, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{X,Y,Tuple{3,default_numlayers(X, Y)}}

                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,nothing},4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y)
                        }
                    }
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{3,nothing,nothing,nothing},4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{3,default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,3,nothing,nothing},4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),3,default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(RightTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,3},4}(
                        unstable_f,
                        Any[-1, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        ccdf,
                        right_fallback
                    )) isa RightTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y),3}
                    }
                end
            end
        end
        @testset "TwoTailCompositeZiggurat" begin
            @testset "Errors" begin
                @test_throws MethodError TwoTailCompositeZiggurat{Float32,Float16}(unstable_f, [-Inf, 0, Inf])
                @test_throws MethodError TwoTailCompositeZiggurat{Float32,Float16}(unstable_f, Any[-Inf, 0, Inf])
                @test_throws MethodError TwoTailCompositeZiggurat{Float32,Float16,nothing}(
                    unstable_f,
                    Any[-Inf, 0, Inf]
                )
                @test_throws MethodError TwoTailCompositeZiggurat{Float32,Float16,3}(unstable_f, Any[-Inf, 0, Inf])
                @test_throws MethodError TwoTailCompositeZiggurat{Float32,Float16,4}(unstable_f, Any[-Inf, 0, Inf])
                @test_throws MethodError TwoTailCompositeZiggurat{Float32,Float16,(nothing, nothing)}(
                    unstable_f,
                    Any[-Inf, 0, Inf]
                )

                @test_throws ErrorException TwoTailCompositeZiggurat{Float32,Float16,Tuple{nothing,nothing},3}(
                    unstable_f,
                    Any[-Inf, 0, Inf]
                )
                @test_throws ErrorException TwoTailCompositeZiggurat{Float32,Float16,Tuple{nothing,nothing},1}(
                    unstable_f,
                    Any[-Inf, 0, Inf]
                )
                @test_throws ErrorException TwoTailCompositeZiggurat{Float32,Float16,NTuple{4,nothing},5}(
                    unstable_f,
                    Any[-Inf, 0, 0.5, 1, Inf]
                )
                @test_throws ErrorException TwoTailCompositeZiggurat{Float32,Float16,NTuple{4,nothing},3}(
                    unstable_f,
                    Any[-Inf, 0, 0.5, 1, Inf]
                )
            end
            @testset for X in (Float64, Float32), Y in (Float64, Float32)
                @testset for (ipdf_left, ipdf_right) in ((nothing, nothing), (unstable_ipdf_left, unstable_ipdf_right)),
                    (cdf, ccdf) in ((nothing, nothing), (unstable_cdf, unstable_ccdf)),
                    (left_fallback, right_fallback) in
                    ((nothing, nothing), (unstable_fallback_left, unstable_fallback_right))

                    # Two parameter forms
                    @test @inferred(TwoTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        (-Inf, 0, Inf);
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,NTuple{2,default_numlayers(X, Y)}}

                    @test @inferred(TwoTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        SVector{3,Any}(-Inf, 0, Inf);
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,NTuple{2,default_numlayers(X, Y)}}

                    @test @inferred(TwoTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        (-Inf, 0, 0.5, Inf);
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,NTuple{3,default_numlayers(X, Y)}}

                    @test @inferred(TwoTailCompositeZiggurat{X,Y}(
                        unstable_f,
                        SVector{4,Any}(-Inf, 0, 0.5, Inf);
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,NTuple{3,default_numlayers(X, Y)}}

                    # Three parameter forms
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{3,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,Tuple{3,default_numlayers(X, Y)}}
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,3}}(
                        unstable_f,
                        Any[-Inf, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),3}}

                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y)
                        }
                    }
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{3,nothing,nothing,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{3,default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,3,nothing,nothing}}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),3,default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,3}}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y),3}
                    }

                    # Four parameter forms
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,nothing,2}(
                        unstable_f,
                        Any[-Inf, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}

                    @test @inferred(TwoTailCompositeZiggurat{X,Y,nothing,4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,NTuple{4,default_numlayers(X, Y)}}

                    @test @inferred(TwoTailCompositeZiggurat{X,Y,3,2}(
                        unstable_f,
                        Any[-Inf, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,Tuple{3,3}}

                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,nothing},2}(
                        unstable_f,
                        Any[-Inf, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,Tuple{default_numlayers(X, Y),default_numlayers(X, Y)}}
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{3,nothing},2}(
                        unstable_f,
                        Any[-Inf, 0, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{X,Y,Tuple{3,default_numlayers(X, Y)}}

                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,nothing},4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y),
                            default_numlayers(X, Y)
                        }
                    }
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{3,nothing,nothing,nothing},4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{3,default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,3,nothing,nothing},4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),3,default_numlayers(X, Y),default_numlayers(X, Y)}
                    }
                    @test @inferred(TwoTailCompositeZiggurat{X,Y,Tuple{nothing,nothing,nothing,3},4}(
                        unstable_f,
                        Any[-Inf, 0, 0.5, 1, Inf];
                        ipdfs = Any[ipdf_left, ipdf_right, ipdf_right, ipdf_right],
                        cdf,
                        ccdf,
                        left_fallback,
                        right_fallback
                    )) isa TwoTailCompositeZiggurat{
                        X,
                        Y,
                        Tuple{default_numlayers(X, Y),default_numlayers(X, Y),default_numlayers(X, Y),3}
                    }
                end
            end
        end
    end
end

nothing
