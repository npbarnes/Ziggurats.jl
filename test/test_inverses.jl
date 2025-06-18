@testset "Generalized Inverses" begin
    @testset "inversepdf()" begin
        import ZigguratTools: inversepdf
        # 0 -> 1
        heaviside1(x) = x >= 0 ? 1.0 : 0.0

        # 0 -> 0
        heaviside2(x) = x > 0 ? 1.0 : 0.0

        # 0 -> 0.5
        function heaviside3(x)
            if x < 0
                0.0
            elseif x > 0
                1.0
            else
                0.5
            end
        end

        function s_curve(x)
            if x < -1
                zero(x)
            elseif x > 1
                2oneunit(x)
            else
                x + oneunit(x)
            end
        end

        # Decreasing functions
        @test inversepdf(cos, (0, π/2))(1/√2) ≈ π / 4

        @test_throws "No inverse exists" inversepdf(cos, (0, π))(2)
        @test_throws "No inverse exists" inversepdf(cos, (0, π))(-2)
        @test inversepdf(cos, (0, π/4))(0) == π/4

        @test inversepdf(x -> s_curve(-x), (-2, 2))(0) >= 1
        @test inversepdf(x -> s_curve(-x), (-2, 2))(1) ≈ 0 atol=1e-16
        @test -2 <= inversepdf(x -> s_curve(-x), (-2, 2))(2) <= -1

        @test inversepdf(x -> heaviside1(-x), (-1, 1))(0) > 0
        @test inversepdf(x -> heaviside1(-x), (-1, 1))(0.5) ≈ 0.0 atol=1e-16
        @test -1 <= inversepdf(x -> heaviside1(-x), (-1, 1))(1) <= 0.0

        @test inversepdf(x -> heaviside2(-x), (-1, 1))(0) >= 0
        @test inversepdf(x -> heaviside2(-x), (-1, 1))(0.5) ≈ 0.0 atol=1e-16
        @test -1 <= inversepdf(x -> heaviside2(-x), (-1, 1))(1) < 0

        @test inversepdf(x -> heaviside3(-x), (-1, 1))(0) > 0
        @test inversepdf(x -> heaviside3(-x), (-1, 1))(0.25) ≈ 0 atol=1e-16
        @test inversepdf(x -> heaviside3(-x), (-1, 1))(0.5) ≈ 0 atol=1e-16
        @test inversepdf(x -> heaviside3(-x), (-1, 1))(0.75) ≈ 0 atol=1e-16
        @test inversepdf(x -> heaviside3(-x), (-1, 1))(1) < 0

        # Increasing functions
        @test inversepdf(cos, (-π/2, 0))(1/√2) ≈ -π / 4

        @test_throws "No inverse exists" inversepdf(cos, (-π, 0))(2)
        @test_throws "No inverse exists" inversepdf(cos, (-π, 0))(-2)
        @test inversepdf(cos, (-π/4, 0))(0) ≈ -π/4

        @test inversepdf(s_curve, (-2, 2))(0) <= -1
        @test inversepdf(s_curve, (-2, 2))(1) ≈ 0 atol=1e-16
        @test inversepdf(s_curve, (-2, 2))(2) >= 1

        @test inversepdf(heaviside1, (-1, 1))(0) < 0
        @test inversepdf(heaviside1, (-1, 1))(0.5) ≈ 0 atol=1e-16
        @test inversepdf(heaviside1, (-1, 1))(1) >= 0

        @test inversepdf(heaviside2, (-1, 1))(0) <= 0
        @test inversepdf(heaviside2, (-1, 1))(0.5) ≈ 0 atol=1e-16
        @test inversepdf(heaviside2, (-1, 1))(1) > 0

        @test inversepdf(heaviside3, (-1, 1))(0) < 0
        @test inversepdf(heaviside3, (-1, 1))(0.25) ≈ 0 atol=1e-16
        @test inversepdf(heaviside3, (-1, 1))(0.5) ≈ 0 atol=1e-16
        @test inversepdf(heaviside3, (-1, 1))(0.75) ≈ 0 atol=1e-16
        @test inversepdf(heaviside3, (-1, 1))(1) > 0

        # Non-monotonicity is not detected
        @test -π/6 <= inversepdf(cos, (-π/6, π/2))(1/√2) <= π/2

        # Constant functions are allowed.
        @test -1 <= inversepdf(x -> 1.0, (-1, 1))(1) <= 1

        # Domains that include some positive numbers
        @test_throws "ArgumentError: No inverse" inversepdf(heaviside3, (-Inf, Inf))(2.0)
        @test_throws "ArgumentError: No inverse" inversepdf(heaviside3, (-Inf, 1.0))(2.0)
        @test_throws "ArgumentError: No inverse" inversepdf(heaviside3, (-1.0, Inf))(2.0)
        @test_throws "ArgumentError: No inverse" inversepdf(heaviside3, (-1.0, 1.0))(2.0)

        @test_throws "ArgumentError: No inverse" inversepdf(heaviside3, (-Inf, Inf))(-2.0)
        @test_throws "ArgumentError: No inverse" inversepdf(heaviside3, (-Inf, 1.0))(-2.0)
        @test_throws "ArgumentError: No inverse" inversepdf(heaviside3, (-1.0, Inf))(-2.0)
        @test_throws "ArgumentError: No inverse" inversepdf(heaviside3, (-1.0, 1.0))(-2.0)

        @test inversepdf(heaviside3, (-Inf, Inf))(1.0) > 0
        @test inversepdf(heaviside3, (-Inf, 1.0))(1.0) > 0
        @test inversepdf(heaviside3, (-1.0, Inf))(1.0) > 0
        @test inversepdf(heaviside3, (-1.0, 1.0))(1.0) > 0

        @test inversepdf(heaviside3, (-Inf, Inf))(0.5) == 0.0
        @test inversepdf(heaviside3, (-Inf, 1.0))(0.5) == 0.0
        @test inversepdf(heaviside3, (-1.0, Inf))(0.5) == 0.0
        @test inversepdf(heaviside3, (-1.0, 1.0))(0.5) == 0.0

        @test inversepdf(heaviside3, (-Inf, Inf))(0) == -Inf
        @test inversepdf(heaviside3, (-Inf, 1.0))(0) == -Inf
        @test inversepdf(heaviside3, (-1.0, Inf))(0) == -1.0
        @test inversepdf(heaviside3, (-1.0, 1.0))(0) == -1.0

        @test_throws "No inverse exists" inversepdf(heaviside3, (-10.0, 10.0))(2.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-10.0, 10.0))(-2.0)
        @test inversepdf(heaviside3, (-10.0, 10.0))(1.0) > 0
        @test inversepdf(heaviside3, (-10.0, 10.0))(0.5) ≈ 0.0 atol = 1e-16
        @test inversepdf(heaviside3, (-10.0, 10.0))(0) < 0

        # domain ends at zero
        @test_throws "No inverse exists" inversepdf(heaviside3, (-Inf, 0.0))(2.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-1.0, 0.0))(2.0)

        @test_throws "No inverse exists" inversepdf(heaviside3, (-Inf, 0.0))(1.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-1.0, 0.0))(1.0)

        @test_throws "No inverse exists" inversepdf(heaviside3, (-Inf, 0.0))(-2.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-1.0, 0.0))(-2.0)

        @test_throws "No inverse exists" inversepdf(heaviside3, (-Inf, 0.0))(-1.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-1.0, 0.0))(-1.0)

        @test inversepdf(heaviside3, (-Inf, 0.0))(0.5) == 0.0
        @test inversepdf(heaviside3, (-1.0, 0.0))(0.5) == 0.0

        @test inversepdf(heaviside3, (-Inf, 0.0))(0) == -Inf
        @test inversepdf(heaviside3, (-1.0, 0.0))(0) == -1.0

        @test_throws "No inverse exists" inversepdf(heaviside3, (-Inf, 0.0))(-2.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-1.0, 0.0))(-2.0)

        @test_throws "No inverse exists" inversepdf(heaviside3, (-10.0, 0.0))(2.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-10.0, 0.0))(1.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-10.0, 0.0))(-2.0)
        @test_throws "No inverse exists" inversepdf(heaviside3, (-10.0, 0.0))(-1.0)
        @test inversepdf(heaviside3, (-10.0, 0.0))(0.5) == 0.0
        @test inversepdf(heaviside3, (-10.0, 0.0))(0) == -10
        @test_throws "No inverse exists" inversepdf(heaviside3, (-10.0, 0.0))(-2.0)

        mheaviside3 = x -> heaviside3(-x)
        @test_throws "No inverse exists" inversepdf(mheaviside3, (-10.0, 10.0))(2.0)
        @test inversepdf(mheaviside3, (-10.0, 10.0))(1.0) < 0.0
        @test inversepdf(mheaviside3, (-10.0, 10.0))(0.75) ≈ 0.0 atol = 1e-16
        @test inversepdf(mheaviside3, (-10.0, 10.0))(0.5) ≈ 0.0 atol = 1e-16
        @test inversepdf(mheaviside3, (-10.0, 10.0))(0.25) ≈ 0.0 atol = 1e-16
        @test inversepdf(mheaviside3, (-10.0, 10.0))(0) > 0.0
        @test_throws "No inverse exists" inversepdf(mheaviside3, (-10.0, 10.0))(-0.25)
        @test_throws "No inverse exists" inversepdf(mheaviside3, (-10.0, 10.0))(-2.0)
    end
end

nothing
