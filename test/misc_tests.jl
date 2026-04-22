@testset "guess_ytype" begin
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-Inf, -12))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-Inf, 0))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-Inf, 12))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-Inf, Inf))) == Float64

    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-20, -12))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-20, 0))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-20, 12))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-20, Inf))) == Float64

    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-5, 0))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-5, 12))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((-5, Inf))) == Float64

    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((5, 8))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((5, 12))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((5, Inf))) == Float64

    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((12, 20))) == Float64
    @test Ziggurats.guess_ytype(exp, Ziggurats.regularize((12, Inf))) == Float64
end
