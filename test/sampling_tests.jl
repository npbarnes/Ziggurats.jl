function monotonic_sampling_tests(td::MonotonicTestData)
    @testset for N in [1, 2, 3, 4]
        @testset "rng isa $(typeof(_rng))" for _rng in [missing, Xoshiro(1234), MersenneTwister(1234)]
            @testset for array_generation in [true, false]
                test_fallbacks = unique([nothing, td.fallback])
                @testset "fallback provided: $(fallback !== nothing)" for fallback in test_fallbacks
                    z = monotonic_ziggurat(td.f, td.domain, N; td.ipdf, td.tailarea, fallback)

                    @test z isa td.constructor
                    @test typeof(rand(z)) == eltype(z) == Ziggurats.Ytype(z) == td.T
                    @test eltype(rand(z, 3)) == td.T

                    test_samples(z, td.dist; td.nbins, td.q, rng = deepcopy(_rng), array_generation)
                end
            end
        end
    end
end

function bell_sampling_tests(td::CompositeTestData)
    @testset for N in [1, 2, 3, 4]
        @testset "rng isa $(typeof(_rng))" for _rng in [missing, Xoshiro(1234), MersenneTwister(1234)]
            @testset for array_generation in [true, false]
                test_leftfallbacks = unique([nothing, td.left_fallback])
                test_rightfallbacks = unique([nothing, td.right_fallback])
                @testset for left_fallback in test_leftfallbacks, right_fallback in test_rightfallbacks
                    @testset "Left side domain" begin
                        z_left = BellZiggurat(
                            td.f,
                            (td.domain[1], td.domain[2]),
                            N;
                            ipdf = td.ipdfs[1],
                            td.cdf,
                            fallback = left_fallback
                        )

                        @test typeof(rand(z_left)) == eltype(z_left) == td.T
                        @test eltype(rand(z_left, 3)) == td.T

                        test_samples(z_left, td.dist; td.nbins, td.q, rng = deepcopy(_rng), array_generation)
                    end

                    @testset "Right side domain" begin
                        z_right = BellZiggurat(
                            td.f,
                            (td.domain[2], td.domain[3]),
                            N;
                            ipdf = td.ipdfs[2],
                            td.ccdf,
                            fallback = right_fallback
                        )

                        @test typeof(rand(z_right)) == eltype(z_right) == td.T
                        @test eltype(rand(z_right, 3)) == td.T

                        test_samples(z_right, td.dist; td.nbins, td.q, rng = deepcopy(_rng), array_generation)
                    end
                end
            end
        end
    end
end

function laplace_sampling_tests(td::CompositeTestData)
    @testset for N in [1, 2, 3, 4]
        @testset "rng isa $(typeof(_rng))" for _rng in [missing, Xoshiro(1234), MersenneTwister(1234)]
            @testset for array_generation in [true, false]
                test_leftfallbacks = unique([nothing, td.left_fallback])
                test_rightfallbacks = unique([nothing, td.right_fallback])
                @testset for left_fallback in test_leftfallbacks, right_fallback in test_rightfallbacks
                    @testset "Left side domain" begin
                        z_left = BellZiggurat(
                            x -> exp(x)/2,
                            (td.domain[1], td.domain[2]),
                            N;
                            ipdf = td.ipdfs[1],
                            td.cdf,
                            fallback = left_fallback
                        )

                        @test typeof(rand(z_left)) == eltype(z_left) == td.T
                        @test eltype(rand(z_left, 3)) == td.T

                        test_samples(z_left, td.dist; td.nbins, td.q, rng = deepcopy(_rng), array_generation)
                    end

                    @testset "Right side domain" begin
                        z_right = BellZiggurat(
                            x -> exp(-x)/2,
                            (td.domain[2], td.domain[3]),
                            N;
                            ipdf = td.ipdfs[2],
                            td.ccdf,
                            fallback = right_fallback
                        )

                        @test typeof(rand(z_right)) == eltype(z_right) == td.T
                        @test eltype(rand(z_right, 3)) == td.T

                        test_samples(z_right, td.dist; td.nbins, td.q, rng = deepcopy(_rng), array_generation)
                    end
                end
            end
        end
    end
end

function composite_sampling_tests(td::CompositeTestData)
    @testset for N in [1, 2, 3, 4]
        @testset "rng isa $(typeof(_rng))" for _rng in [missing, Xoshiro(1234), MersenneTwister(1234)]
            @testset for array_generation in [true, false]
                test_leftfallbacks = unique([nothing, td.left_fallback])
                test_rightfallbacks = unique([nothing, td.right_fallback])
                @testset for left_fallback in test_leftfallbacks, right_fallback in test_rightfallbacks
                    z = CompositeZiggurat(td.f, td.domain, N; td.ipdfs, td.cdf, td.ccdf, left_fallback, right_fallback)

                    @test typeof(rand(z)) == eltype(z) == td.T
                    @test eltype(rand(z, 3)) == td.T

                    test_samples(z, td.dist; td.nbins, td.q, rng = deepcopy(_rng), array_generation)
                end
            end
        end
    end
end

@testset "Sampling Tests" begin
    @testset "Monotonic" begin
        @testset "$(td.name)" for td in MonotonicTestCases
            monotonic_sampling_tests(td)
        end
    end

    @testset "Bell" begin
        @testset "$(td.name)" for td in BellTestCases
            bell_sampling_tests(td)
        end
        @testset "Laplace w/ half-pdf" for td in LaplaceTestCases
            laplace_sampling_tests(td)
        end
    end

    @testset "Composite" begin
        @testset "$(td.name)" for td in CompositeTestCases
            composite_sampling_tests(td)
        end
    end
end

nothing
