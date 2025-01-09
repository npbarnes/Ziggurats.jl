@testset "Generalized Inverses" begin
    @testset "_middle()" begin
        import ZigguratTools: inverse, _middle

        # Opposite signs have zero in the middle
        @test _middle(0.0, -0.0) == 0.0
        @test _middle(-1.0, 1.0) == 0.0
        @test _middle(-1.0, nextfloat(0.0)) == 0.0
        @test _middle(prevfloat(0.0), 1.0) == 0.0

        @test _middle(-1.0, Inf) == 0.0
        @test _middle(-Inf, 1.0) == 0.0
        @test _middle(-Inf, Inf) == 0.0
        @test _middle(-Inf, nextfloat(0.0)) == 0.0
        @test _middle(prevfloat(0.0), Inf) == 0.0

        @test _middle(-1.0, prevfloat(Inf)) == 0.0
        @test _middle(nextfloat(-Inf), 1.0) == 0.0
        @test _middle(nextfloat(-Inf), prevfloat(Inf)) == 0.0
        @test _middle(-Inf, prevfloat(Inf)) == 0.0
        @test _middle(nextfloat(-Inf), Inf) == 0.0
        @test _middle(nextfloat(-Inf), nextfloat(0.0)) == 0.0
        @test _middle(prevfloat(0.0), prevfloat(Inf)) == 0.0

        # Edge cases
        @test _middle(prevfloat(prevfloat(Inf)), Inf) == prevfloat(Inf)
        @test _middle(-Inf, nextfloat(nextfloat(-Inf))) == nextfloat(-Inf)
        @test 0 < _middle(-0.0, 1.0) < 1 # negative zero isn't negative

        # Typical cases
        @test 1e100 < _middle(1e100, 1e200) < 1e200
        @test -1e200 < _middle(-1e200, -1e100) < -1e100
        @test _middle(1e200, nextfloat(nextfloat(1e200))) == nextfloat(1e200)
        @test _middle(prevfloat(prevfloat(-1e200)), -1e200) == prevfloat(-1e200)

        # Commutativity
        @test _middle(0.0, -0.0) == _middle(-0.0, 0.0)
        @test _middle(-1.0, 1.0) == _middle(1.0, -1.0)
        @test _middle(-1.0, nextfloat(0.0)) == _middle(nextfloat(0.0), -1.0)
        @test _middle(prevfloat(0.0), 1.0) == _middle(1.0, prevfloat(0.0))

        @test _middle(-1.0, Inf) == _middle(Inf, -1.0)
        @test _middle(-Inf, 1.0) == _middle(1.0, -Inf)
        @test _middle(-Inf, Inf) == _middle(Inf, -Inf)
        @test _middle(-Inf, nextfloat(0.0)) == _middle(nextfloat(0.0), -Inf)
        @test _middle(prevfloat(0.0), Inf) == _middle(Inf, prevfloat(0.0))

        @test _middle(-1.0, prevfloat(Inf)) == _middle(prevfloat(Inf), -1.0)
        @test _middle(nextfloat(-Inf), 1.0) == _middle(1.0, nextfloat(-Inf))
        @test _middle(nextfloat(-Inf), prevfloat(Inf)) ==
              _middle(prevfloat(Inf), nextfloat(-Inf))
        @test _middle(-Inf, prevfloat(Inf)) == _middle(prevfloat(Inf), -Inf)
        @test _middle(nextfloat(-Inf), Inf) == _middle(Inf, nextfloat(-Inf))
        @test _middle(nextfloat(-Inf), nextfloat(0.0)) ==
              _middle(nextfloat(0.0), nextfloat(-Inf))
        @test _middle(prevfloat(0.0), prevfloat(Inf)) ==
              _middle(prevfloat(Inf), prevfloat(0.0))

        @test _middle(prevfloat(prevfloat(Inf)), Inf) ==
              _middle(Inf, prevfloat(prevfloat(Inf)))
        @test _middle(-Inf, nextfloat(nextfloat(-Inf))) ==
              _middle(nextfloat(nextfloat(-Inf)), -Inf)
        @test _middle(-0.0, 1.0) == _middle(1.0, -0.0)

        @test _middle(1e100, 1e200) == _middle(1e200, 1e100)
        @test _middle(-1e200, -1e100) == _middle(-1e100, -1e200)
        @test _middle(1e200, nextfloat(nextfloat(1e200))) ==
              _middle(nextfloat(nextfloat(1e200)), 1e200)
        @test _middle(prevfloat(prevfloat(-1e200)), -1e200) ==
              _middle(-1e200, prevfloat(prevfloat(-1e200)))
    end

    @testset "inverse()" begin
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

        # This function is special because slowdecay(x) > nextfloat(0.0) for all finite x.
        slowdecay(x) = 1 / log(1 + x)

        # Decreasing functions
        @test inverse(cos, (0, π), 0) == π / 2
        @test inverse(slowdecay, (0, Inf), 0) == Inf
        @test inverse(slowdecay, (0, Inf), nextfloat(0.0)) == prevfloat(Inf)
        @test begin
            x = inverse(slowdecay, (0, Inf), Inf)
            if isinf(slowdecay(x))
                !isinf(slowdecay(nextfloat(x)))
            end
        end

        @test_throws "no solutions exist" inverse(cos, (0, π), 2)
        @test inverse(cos, (0, π), -2) == float(π)
        @test inverse(slowdecay, (0, Inf), -2) == Inf
        @test inverse(slowdecay, (0, Inf), nextfloat(0.0)) == prevfloat(Inf)

        @test inverse(x -> s_curve(-x), (-2, 2), -1) == 2
        @test inverse(x -> s_curve(-x), (-2, 2), 0) == 0
        @test inverse(x -> s_curve(-x), (-2, 2), 1) == -1

        @test inverse(x -> heaviside1(-x), (-1, 1), 0) == 1
        @test inverse(x -> heaviside1(-x), (-1, 1), 0.5) == 0.0
        @test inverse(x -> heaviside1(-x), (-1, 1), 1) == 0.0

        @test inverse(x -> heaviside2(-x), (-1, 1), 0) == 1
        @test inverse(x -> heaviside2(-x), (-1, 1), 0.5) == prevfloat(0.0)
        @test inverse(x -> heaviside2(-x), (-1, 1), 1) == prevfloat(0.0)

        @test inverse(x -> heaviside3(-x), (-1, 1), 0) == 1
        @test inverse(x -> heaviside3(-x), (-1, 1), 0.25) == 0.0
        @test inverse(x -> heaviside3(-x), (-1, 1), 0.5) == 0.0
        @test inverse(x -> heaviside3(-x), (-1, 1), 0.75) == prevfloat(0.0)
        @test inverse(x -> heaviside3(-x), (-1, 1), 1) == prevfloat(0.0)

        # Increasing functions
        @test inverse(cos, (-π, 0), 0) == -π / 2
        @test inverse(x -> slowdecay(-x), (-Inf, 0), 0) == -Inf
        @test inverse(x -> slowdecay(-x), (-Inf, 0), nextfloat(0.0)) == nextfloat(-Inf)

        @test_throws "no solutions exist" inverse(cos, (-π, 0), 2)
        @test inverse(cos, (-π, 0), -2) == float(-π)

        @test inverse(s_curve, (-2, 2), -1) == -2
        @test inverse(s_curve, (-2, 2), 0) == 0
        @test inverse(s_curve, (-2, 2), 1) == 1

        @test inverse(heaviside1, (-1, 1), 0) == -1
        @test inverse(heaviside1, (-1, 1), 0.5) == 0.0
        @test inverse(heaviside1, (-1, 1), 1) == 0.0

        @test inverse(heaviside2, (-1, 1), 0) == -1
        @test inverse(heaviside2, (-1, 1), 0.5) == nextfloat(0.0)
        @test inverse(heaviside2, (-1, 1), 1) == nextfloat(0.0)

        @test inverse(heaviside3, (-1, 1), 0) == -1
        @test inverse(heaviside3, (-1, 1), 0.25) == 0.0
        @test inverse(heaviside3, (-1, 1), 0.5) == 0.0
        @test inverse(heaviside3, (-1, 1), 0.75) == nextfloat(0.0)
        @test inverse(heaviside3, (-1, 1), 1) == nextfloat(0.0)

        # Non-monotonicity is sometimes detected
        @test_throws "f must be monotonic" inverse(cos, (-1.3, π), 0)

        # Non-monotonicity is not always detected. Ideally, it would be, but that's
        # probably impossible to do in general.
        @test_broken try
            inverse(x -> x == 1 ? 7.0 : cos(x), (0, π), 0)
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
        @test_throws "f must be non-constant" inverse(x -> 1.0, (-1, 1), 1)
        @test_throws "f must be non-constant" inverse(sign, (-2, -1), -1)
        @test_throws "f must be non-constant" inverse(sign, (1, 2), 1)

        # These have no solution, but f being constant is detected first.
        @test_throws "f must be non-constant" inverse(x -> 1.0, (-1, 1), 2)
        @test_throws "f must be non-constant" inverse(sign, (-2, -1), -2)
        @test_throws "f must be non-constant" inverse(sign, (1, 2), 2)

        # Sometimes constant functions and non-monotonic functions get mixed up.
        @test_broken try
            inverse(cos, (-1, 1), 0)
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

        @test_throws "domain must be an ordered tuple" inverse(cos, (π, 0), 0)
        @test_throws "domain must be an ordered tuple" inverse(s_curve, (1, -1), 0)
        @test_throws "domain must be an ordered tuple" inverse(slowdecay, (Inf, 0), 0)

        # Domains that include some positive numbers
        @test_throws "no solutions" inverse(sign, (-Inf, Inf), 2.0)
        @test_throws "no solutions" inverse(sign, (-Inf, 1.0), 2.0)
        @test_throws "no solutions" inverse(sign, (-1.0, Inf), 2.0)
        @test_throws "no solutions" inverse(sign, (-1.0, 1.0), 2.0)

        @test inverse(sign, (-Inf, Inf), 1.0) == nextfloat(0.0)
        @test inverse(sign, (-Inf, 1.0), 1.0) == nextfloat(0.0)
        @test inverse(sign, (-1.0, Inf), 1.0) == nextfloat(0.0)
        @test inverse(sign, (-1.0, 1.0), 1.0) == nextfloat(0.0)

        @test inverse(sign, (-Inf, Inf), 0.0) == 0.0
        @test inverse(sign, (-Inf, 1.0), 0.0) == 0.0
        @test inverse(sign, (-1.0, Inf), 0.0) == 0.0
        @test inverse(sign, (-1.0, 1.0), 0.0) == 0.0

        @test inverse(sign, (-Inf, Inf), -1.0) == -Inf
        @test inverse(sign, (-Inf, 1.0), -1.0) == -Inf
        @test inverse(sign, (-1.0, Inf), -1.0) == -1.0
        @test inverse(sign, (-1.0, 1.0), -1.0) == -1.0

        @test inverse(sign, (-Inf, Inf), -2.0) == -Inf
        @test inverse(sign, (-Inf, 1.0), -2.0) == -Inf
        @test inverse(sign, (-1.0, Inf), -2.0) == -1.0
        @test inverse(sign, (-1.0, 1.0), -2.0) == -1.0

        @test_throws "no solutions exist" inverse(sign, (-10.0, 10.0), 2.0)
        @test inverse(sign, (-10.0, 10.0), 1.0) == nextfloat(0.0)
        @test inverse(sign, (-10.0, 10.0), 0.0) == 0.0
        @test inverse(sign, (-10.0, 10.0), -1.0) == -10
        @test inverse(sign, (-10.0, 10.0), -2.0) == -10.0

        # domain ends at zero
        @test_throws "no solutions exist" inverse(sign, (-Inf, 0.0), 2.0)
        @test_throws "no solutions exist" inverse(sign, (-1.0, 0.0), 2.0)

        @test_throws "no solutions exist" inverse(sign, (-Inf, 0.0), 1.0)
        @test_throws "no solutions exist" inverse(sign, (-1.0, 0.0), 1.0)

        @test inverse(sign, (-Inf, 0.0), 0.0) == 0.0
        @test inverse(sign, (-1.0, 0.0), 0.0) == 0.0

        @test inverse(sign, (-Inf, 0.0), -1.0) == -Inf
        @test inverse(sign, (-1.0, 0.0), -1.0) == -1.0

        @test inverse(sign, (-Inf, 0.0), -2.0) == -Inf
        @test inverse(sign, (-1.0, 0.0), -2.0) == -1.0

        @test_throws "no solutions exist" inverse(sign, (-10.0, 0.0), 2.0)
        @test_throws "no solutions exist" inverse(sign, (-10.0, 0.0), 1.0)
        @test inverse(sign, (-10.0, 0.0), 0.0) == 0.0
        @test inverse(sign, (-10.0, 0.0), -1.0) == -10
        @test inverse(sign, (-10.0, 0.0), -2.0) == -10.0

        # Out of order domains
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -1.0), 2.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -1.0), 2.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -Inf), 2.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -Inf), 2.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -1.0), 1.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -1.0), 1.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -Inf), 1.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -Inf), 1.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -1.0), 0.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -1.0), 0.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -Inf), 0.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -Inf), 0.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -1.0), -1.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -1.0), -1.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -Inf), -1.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -Inf), -1.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -1.0), -2.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -1.0), -2.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (1.0, -Inf), -2.0)
        @test_throws "domain must be an ordered tuple" inverse(sign, (Inf, -Inf), -2.0)

        msign = x -> -sign(x)
        @test_throws "no solutions exist" inverse(msign, (-10.0, 10.0), 2.0)
        @test inverse(msign, (-10.0, 10.0), 1.0) == prevfloat(0.0)
        @test inverse(msign, (-10.0, 10.0), 0.5) == prevfloat(0.0)
        @test inverse(msign, (-10.0, 10.0), 0.0) == 0.0
        @test inverse(msign, (-10.0, 10.0), -0.5) == 0.0
        @test inverse(msign, (-10.0, 10.0), -1.0) == 10.0
        @test inverse(msign, (-10.0, 10.0), -2.0) == 10.0
    end

    @testset "Handle Floating Point Non-Monotone with xatol" begin
        # Chisq(3) distribution in the form exp(logpdf(Chisq(3)))
        # The extra set of paretheses is significant.
        using SpecialFunctions
        f(x) = exp((-log(gamma(1.5)) - log(2) - x / 2) + (1.5 - 1) * log(x / 2))
        x = 0.19761117965603436
        y = 0.16065901454232936

        # This function is mathematically monotonic on (0,1), but in floating
        # point this implementation is not monotonic.
        @test_throws "f must be monotonic" inverse(f, (0, 1), y)
        @test f(inverse(f, (0, 1), y; xatol = eps(1.0))) ≈ y
    end
end

nothing
