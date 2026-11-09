@testset "guess_ytype" begin
    @test Ziggurats.guess_ytype(exp, (-Inf, -12.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (-Inf, 0.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (-Inf, 12.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (-Inf, Inf)) == Float64

    @test Ziggurats.guess_ytype(exp, (-20.0, -12.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (-20.0, 0.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (-20.0, 12.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (-20.0, Inf)) == Float64

    @test Ziggurats.guess_ytype(exp, (-5.0, 0.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (-5.0, 12.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (-5.0, Inf)) == Float64

    @test Ziggurats.guess_ytype(exp, (5.0, 8.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (5.0, 12.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (5.0, Inf)) == Float64

    @test Ziggurats.guess_ytype(exp, (12.0, 20.0)) == Float64
    @test Ziggurats.guess_ytype(exp, (12.0, Inf)) == Float64
end
