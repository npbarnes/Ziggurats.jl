issamevector(a::Vector{T}, b::Vector{T}) where {T} = a == b
issamevector(a::Vector, b::Vector) = false

@testset "Argument handling" begin
    @testset "regularize domains" begin
        import ZigguratTools: regularize
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
end
