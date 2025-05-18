@testset "Generalized Inverses" begin
    @testset "inverse()" begin
        import ZigguratTools: inverse
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
                -oneunit(x)
            elseif x > 1
                oneunit(x)
            else
                x
            end
        end

        # Decreasing functions
        @test inverse(cos, (0, π))(0) == π / 2

        @test_throws "No inverse exists" inverse(cos, (0, π))(2)
        @test inverse(cos, (0, π))(-2) == float(π)

        @test inverse(x -> s_curve(-x), (-2, 2))(-1) == 2
        @test inverse(x -> s_curve(-x), (-2, 2))(0) == 0
        @test inverse(x -> s_curve(-x), (-2, 2))(1) == -1

        @test inverse(x -> heaviside1(-x), (-1, 1))(0) == 1
        @test inverse(x -> heaviside1(-x), (-1, 1))(0.5) == 0.0
        @test inverse(x -> heaviside1(-x), (-1, 1))(1) == 0.0

        @test inverse(x -> heaviside2(-x), (-1, 1))(0) == 1
        @test inverse(x -> heaviside2(-x), (-1, 1))(0.5) == prevfloat(0.0)
        @test inverse(x -> heaviside2(-x), (-1, 1))(1) == prevfloat(0.0)

        @test inverse(x -> heaviside3(-x), (-1, 1))(0) == 1
        @test inverse(x -> heaviside3(-x), (-1, 1))(0.25) == 0.0
        @test inverse(x -> heaviside3(-x), (-1, 1))(0.5) == 0.0
        @test inverse(x -> heaviside3(-x), (-1, 1))(0.75) == prevfloat(0.0)
        @test inverse(x -> heaviside3(-x), (-1, 1))(1) == prevfloat(0.0)

        # Increasing functions
        @test inverse(cos, (-π, 0))(0) == -π / 2

        @test_throws "No inverse exists" inverse(cos, (-π, 0))(2)
        @test inverse(cos, (-π, 0))(-2) == float(-π)

        @test inverse(s_curve, (-2, 2))(-1) == -2
        @test inverse(s_curve, (-2, 2))(0) == 0
        @test inverse(s_curve, (-2, 2))(1) == 1

        @test inverse(heaviside1, (-1, 1))(0) == -1
        @test inverse(heaviside1, (-1, 1))(0.5) == 0.0
        @test inverse(heaviside1, (-1, 1))(1) == 0.0

        @test inverse(heaviside2, (-1, 1))(0) == -1
        @test inverse(heaviside2, (-1, 1))(0.5) == nextfloat(0.0)
        @test inverse(heaviside2, (-1, 1))(1) == nextfloat(0.0)

        @test inverse(heaviside3, (-1, 1))(0) == -1
        @test inverse(heaviside3, (-1, 1))(0.25) == 0.0
        @test inverse(heaviside3, (-1, 1))(0.5) == 0.0
        @test inverse(heaviside3, (-1, 1))(0.75) == nextfloat(0.0)
        @test inverse(heaviside3, (-1, 1))(1) == nextfloat(0.0)

        # Non-monotonicity is sometimes detected
        @test_throws "f must be monotonic" inverse(cos, (-1.3, π))(0)

        # Non-monotonicity is not always detected. Ideally, it would be, but that's
        # probably impossible to do in general.
        @test_broken try
            inverse(x -> x == 1 ? 7.0 : cos(x), (0, π))(0)
        catch e
            if e.msg == "f must be monotonic"
                true # test passes if non-monotonicity is detected.
            else
                rethrow()
            end
        else
            false # test fails if no error is thrown
        end

        # Constant functions are not allowed.
        @test_throws "f must be non-constant" inverse(x -> 1.0, (-1, 1))(1)
        @test_throws "f must be non-constant" inverse(sign, (-2, -1))(-1)
        @test_throws "f must be non-constant" inverse(sign, (1, 2))(1)

        # These have no solution, but f being constant is detected first.
        @test_throws "f must be non-constant" inverse(x -> 1.0, (-1, 1))(2)
        @test_throws "f must be non-constant" inverse(sign, (-2, -1))(-2)
        @test_throws "f must be non-constant" inverse(sign, (1, 2))(2)

        # Sometimes constant functions and non-monotonic functions get mixed up.
        @test_broken try
            inverse(cos, (-1, 1))(0)
        catch e
            if e.msg == "f must be monotonic"
                true # test passes
            elseif e.msg == "f must be non-constant"
                false # test fails
            else
                rethrow()
            end
        else
            false # test fails if no error is thrown
        end

        # Domains that include some positive numbers
        @test_throws "no solutions" inverse(sign, (-Inf, Inf))(2.0)
        @test_throws "no solutions" inverse(sign, (-Inf, 1.0))(2.0)
        @test_throws "no solutions" inverse(sign, (-1.0, Inf))(2.0)
        @test_throws "no solutions" inverse(sign, (-1.0, 1.0))(2.0)

        @test inverse(sign, (-Inf, Inf))(1.0) == nextfloat(0.0)
        @test inverse(sign, (-Inf, 1.0))(1.0) == nextfloat(0.0)
        @test inverse(sign, (-1.0, Inf))(1.0) == nextfloat(0.0)
        @test inverse(sign, (-1.0, 1.0))(1.0) == nextfloat(0.0)

        @test inverse(sign, (-Inf, Inf))(0.0) == 0.0
        @test inverse(sign, (-Inf, 1.0))(0.0) == 0.0
        @test inverse(sign, (-1.0, Inf))(0.0) == 0.0
        @test inverse(sign, (-1.0, 1.0))(0.0) == 0.0

        @test inverse(sign, (-Inf, Inf))(-1.0) == -Inf
        @test inverse(sign, (-Inf, 1.0))(-1.0) == -Inf
        @test inverse(sign, (-1.0, Inf))(-1.0) == -1.0
        @test inverse(sign, (-1.0, 1.0))(-1.0) == -1.0

        @test inverse(sign, (-Inf, Inf))(-2.0) == -Inf
        @test inverse(sign, (-Inf, 1.0))(-2.0) == -Inf
        @test inverse(sign, (-1.0, Inf))(-2.0) == -1.0
        @test inverse(sign, (-1.0, 1.0))(-2.0) == -1.0

        @test_throws "No inverse exists" inverse(sign, (-10.0, 10.0))(2.0)
        @test inverse(sign, (-10.0, 10.0))(1.0) == nextfloat(0.0)
        @test inverse(sign, (-10.0, 10.0))(0.0) == 0.0
        @test inverse(sign, (-10.0, 10.0))(-1.0) == -10
        @test inverse(sign, (-10.0, 10.0))(-2.0) == -10.0

        # domain ends at zero
        @test_throws "No inverse exists" inverse(sign, (-Inf, 0.0))(2.0)
        @test_throws "No inverse exists" inverse(sign, (-1.0, 0.0))(2.0)

        @test_throws "No inverse exists" inverse(sign, (-Inf, 0.0))(1.0)
        @test_throws "No inverse exists" inverse(sign, (-1.0, 0.0))(1.0)

        @test inverse(sign, (-Inf, 0.0))(0.0) == 0.0
        @test inverse(sign, (-1.0, 0.0))(0.0) == 0.0

        @test inverse(sign, (-Inf, 0.0))(-1.0) == -Inf
        @test inverse(sign, (-1.0, 0.0))(-1.0) == -1.0

        @test inverse(sign, (-Inf, 0.0))(-2.0) == -Inf
        @test inverse(sign, (-1.0, 0.0))(-2.0) == -1.0

        @test_throws "No inverse exists" inverse(sign, (-10.0, 0.0))(2.0)
        @test_throws "No inverse exists" inverse(sign, (-10.0, 0.0))(1.0)
        @test inverse(sign, (-10.0, 0.0))(0.0) == 0.0
        @test inverse(sign, (-10.0, 0.0))(-1.0) == -10
        @test inverse(sign, (-10.0, 0.0))(-2.0) == -10.0

        msign = x -> -sign(x)
        @test_throws "No inverse exists" inverse(msign, (-10.0, 10.0))(2.0)
        @test inverse(msign, (-10.0, 10.0))(1.0) == prevfloat(0.0)
        @test inverse(msign, (-10.0, 10.0))(0.5) == prevfloat(0.0)
        @test inverse(msign, (-10.0, 10.0))(0.0) == 0.0
        @test inverse(msign, (-10.0, 10.0))(-0.5) == 0.0
        @test inverse(msign, (-10.0, 10.0))(-1.0) == 10.0
        @test inverse(msign, (-10.0, 10.0))(-2.0) == 10.0
    end
end

nothing
