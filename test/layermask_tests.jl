@testset "layermask" begin
    struct Foo end
    struct Bar <: Number end
    @testset "without sign bit" begin
        import Ziggurats: layermask
        @testset "Error when the number of layers is non-positive" begin
            @test_throws ArgumentError layermask(Float16, 0)
            @test_throws ArgumentError layermask(Float16, -1)
            @test_throws ArgumentError layermask(Float32, 0)
            @test_throws ArgumentError layermask(Float32, -1)
            @test_throws ArgumentError layermask(Float64, 0)
            @test_throws ArgumentError layermask(Float64, -1)

            # Non-FloatXX also throws
            @test_throws ArgumentError layermask(UInt64, 0)
            @test_throws ArgumentError layermask(UInt64, -1)
            @test_throws ArgumentError layermask(Any, 0)
            @test_throws ArgumentError layermask(Any, -1)
            @test_throws ArgumentError layermask(Number, 0)
            @test_throws ArgumentError layermask(Number, -1)
            @test_throws ArgumentError layermask(Foo, 0)
            @test_throws ArgumentError layermask(Foo, -1)
            @test_throws ArgumentError layermask(Bar, 0)
            @test_throws ArgumentError layermask(Bar, -1)
        end
        @testset "return (nothing, false) when T is not FloatXX and N is positive" begin
            # Including powers of two and non-powers of two and too-large numbers
            @test layermask(UInt64, 1) === (nothing, false)
            @test layermask(UInt64, 10) === (nothing, false)
            @test layermask(UInt64, 256) === (nothing, false)
            @test layermask(UInt64, big"2"^65) === (nothing, false)
            @test layermask(UInt64, big"2"^65 + 1) === (nothing, false)
            @test layermask(Any, 1) === (nothing, false)
            @test layermask(Any, 10) === (nothing, false)
            @test layermask(Any, 256) === (nothing, false)
            @test layermask(Any, big"2"^65) === (nothing, false)
            @test layermask(Any, big"2"^65 + 1) === (nothing, false)
            @test layermask(Number, 1) === (nothing, false)
            @test layermask(Number, 10) === (nothing, false)
            @test layermask(Number, 256) === (nothing, false)
            @test layermask(Number, big"2"^65) === (nothing, false)
            @test layermask(Number, big"2"^65 + 1) === (nothing, false)
            @test layermask(Foo, 1) === (nothing, false)
            @test layermask(Foo, 10) === (nothing, false)
            @test layermask(Foo, 256) === (nothing, false)
            @test layermask(Foo, big"2"^65) === (nothing, false)
            @test layermask(Foo, big"2"^65 + 1) === (nothing, false)
            @test layermask(Bar, 1) === (nothing, false)
            @test layermask(Bar, 10) === (nothing, false)
            @test layermask(Bar, 256) === (nothing, false)
            @test layermask(Bar, big"2"^65) === (nothing, false)
            @test layermask(Bar, big"2"^65 + 1) === (nothing, false)
        end
        @testset "return (nothing, false) when N is too large for a bitmask" begin
            # Powers of two
            @test layermask(Float64, big(2)^65) === (nothing, false)
            @test layermask(Float32, 2^33) === (nothing, false)
            @test layermask(Float16, 2^17) === (nothing, false)

            # Non-powers of two
            @test layermask(Float64, big(2)^64 + 1) === (nothing, false)
            @test layermask(Float32, 2^32 + 1) === (nothing, false)
            @test layermask(Float16, 2^16 + 1) === (nothing, false)
        end
        @testset "return (nothing, false) when N is not a power of 2" begin
            @test layermask(Float64, 100) === (nothing, false)
            @test layermask(Float32, 100) === (nothing, false)
            @test layermask(Float16, 100) === (nothing, false)
        end
        @testset "when N = 2^m the bitmask has the first m bits set" begin
            # when N = 2^m is small enough and T isa FloatXX, then the result is a bitmask with
            # the same bitwidth as T with the first m least significant bits set.

            # Float16 → UInt16
            @test layermask(Float16, 1)[1] === UInt16(0x0000)
            @test layermask(Float16, 2)[1] === UInt16(0x0001)
            @test layermask(Float16, 4)[1] === UInt16(0x0003)
            @test layermask(Float16, 8)[1] === UInt16(0x0007)
            @test layermask(Float16, 2^16)[1] === UInt16(0xFFFF)

            # Float32 → UInt32
            @test layermask(Float32, 1)[1] === UInt32(0x00000000)
            @test layermask(Float32, 2)[1] === UInt32(0x00000001)
            @test layermask(Float32, 4)[1] === UInt32(0x00000003)
            @test layermask(Float32, 8)[1] === UInt32(0x00000007)
            @test layermask(Float32, 2^32)[1] === UInt32(0xFFFFFFFF)

            # Float64 → UInt64
            @test layermask(Float64, 1)[1] === UInt64(0x0000000000000000)
            @test layermask(Float64, 2)[1] === UInt64(0x0000000000000001)
            @test layermask(Float64, 4)[1] === UInt64(0x0000000000000003)
            @test layermask(Float64, 8)[1] === UInt64(0x0000000000000007)
            @test layermask(Float64, big(2)^64)[1] === UInt64(0xFFFFFFFFFFFFFFFF)
        end
        @testset "second return value indicates overlap" begin
            # Float16 → UInt16
            @test layermask(Float16, 1)[2] === false
            @test layermask(Float16, 2)[2] === false
            @test layermask(Float16, 4)[2] === false
            @test layermask(Float16, 8)[2] === false
            @test layermask(Float16, 2^6)[2] === false
            @test layermask(Float16, 2^7)[2] === true
            @test layermask(Float16, 2^16)[2] === true

            # Float32 → UInt32
            @test layermask(Float32, 1)[2] === false
            @test layermask(Float32, 2)[2] === false
            @test layermask(Float32, 4)[2] === false
            @test layermask(Float32, 8)[2] === false
            @test layermask(Float32, 2^9)[2] === false
            @test layermask(Float32, 2^10)[2] === true
            @test layermask(Float32, 2^32)[2] === true

            # Float64 → UInt64
            @test layermask(Float64, 1)[2] === false
            @test layermask(Float64, 2)[2] === false
            @test layermask(Float64, 4)[2] === false
            @test layermask(Float64, 8)[2] === false
            @test layermask(Float64, 2^12)[2] === false
            @test layermask(Float64, 2^13)[2] === true
            @test layermask(Float64, big(2)^64)[2] === true
        end
    end

    @testset "with sign bit" begin
        import Ziggurats: layermask_signed
        @testset "Error when the number of layers is non-positive" begin
            @test_throws ArgumentError layermask_signed(Float16, 0)
            @test_throws ArgumentError layermask_signed(Float16, -1)
            @test_throws ArgumentError layermask_signed(Float32, 0)
            @test_throws ArgumentError layermask_signed(Float32, -1)
            @test_throws ArgumentError layermask_signed(Float64, 0)
            @test_throws ArgumentError layermask_signed(Float64, -1)

            # Non-FloatXX also throws
            @test_throws ArgumentError layermask_signed(UInt64, 0)
            @test_throws ArgumentError layermask_signed(UInt64, -1)
            @test_throws ArgumentError layermask_signed(Any, 0)
            @test_throws ArgumentError layermask_signed(Any, -1)
            @test_throws ArgumentError layermask_signed(Number, 0)
            @test_throws ArgumentError layermask_signed(Number, -1)
            @test_throws ArgumentError layermask_signed(Foo, 0)
            @test_throws ArgumentError layermask_signed(Foo, -1)
            @test_throws ArgumentError layermask_signed(Bar, 0)
            @test_throws ArgumentError layermask_signed(Bar, -1)
        end
        @testset "return (nothing, false) when T is not FloatXX and N is positive" begin
            # Including powers of two and non-powers of two and too-large numbers
            @test layermask_signed(UInt64, 1) === (nothing, false)
            @test layermask_signed(UInt64, 10) === (nothing, false)
            @test layermask_signed(UInt64, 256) === (nothing, false)
            @test layermask_signed(UInt64, big"2"^65) === (nothing, false)
            @test layermask_signed(UInt64, big"2"^65 + 1) === (nothing, false)
            @test layermask_signed(Any, 1) === (nothing, false)
            @test layermask_signed(Any, 10) === (nothing, false)
            @test layermask_signed(Any, 256) === (nothing, false)
            @test layermask_signed(Any, big"2"^65) === (nothing, false)
            @test layermask_signed(Any, big"2"^65 + 1) === (nothing, false)
            @test layermask_signed(Number, 1) === (nothing, false)
            @test layermask_signed(Number, 10) === (nothing, false)
            @test layermask_signed(Number, 256) === (nothing, false)
            @test layermask_signed(Number, big"2"^65) === (nothing, false)
            @test layermask_signed(Number, big"2"^65 + 1) === (nothing, false)
            @test layermask_signed(Foo, 1) === (nothing, false)
            @test layermask_signed(Foo, 10) === (nothing, false)
            @test layermask_signed(Foo, 256) === (nothing, false)
            @test layermask_signed(Foo, big"2"^65) === (nothing, false)
            @test layermask_signed(Foo, big"2"^65 + 1) === (nothing, false)
            
            @test layermask_signed(Bar, 1) === (nothing, false)
            @test layermask_signed(Bar, 10) === (nothing, false)
            @test layermask_signed(Bar, 256) === (nothing, false)
            @test layermask_signed(Bar, big"2"^65) === (nothing, false)
            @test layermask_signed(Bar, big"2"^65 + 1) === (nothing, false)
        end
        @testset "return (nothing, false) when N is too large for a bitmask" begin
            # Powers of two
            @test layermask_signed(Float64, big(2)^64) === (nothing, false)
            @test layermask_signed(Float32, 2^32) === (nothing, false)
            @test layermask_signed(Float16, 2^16) === (nothing, false)

            # Non-powers of two
            @test layermask_signed(Float64, big(2)^63 + 1) === (nothing, false)
            @test layermask_signed(Float32, 2^31 + 1) === (nothing, false)
            @test layermask_signed(Float16, 2^15 + 1) === (nothing, false)
        end
        @testset "return (nothing, false) when N is not a power of 2" begin
            @test layermask_signed(Float64, 100) === (nothing, false)
            @test layermask_signed(Float32, 100) === (nothing, false)
            @test layermask_signed(Float16, 100) === (nothing, false)
        end
        @testset "when N = 2^m the bitmask has the first m bits set shifted by one" begin
            # when N = 2^m is small enough and T isa FloatXX, then the result is a bitmask with
            # the same bitwidth as T with m consecutive bits set starting with the second bit.

            # Float16 → UInt16
            @test layermask_signed(Float16, 1)[1] === UInt16(0x0000)
            @test layermask_signed(Float16, 2)[1] === UInt16(0x0002)
            @test layermask_signed(Float16, 4)[1] === UInt16(0x0006)
            @test layermask_signed(Float16, 8)[1] === UInt16(0x000E)
            @test layermask_signed(Float16, 2^15)[1] === UInt16(0xFFFE)

            # Float32 → UInt32
            @test layermask_signed(Float32, 1)[1] === UInt32(0x00000000)
            @test layermask_signed(Float32, 2)[1] === UInt32(0x00000002)
            @test layermask_signed(Float32, 4)[1] === UInt32(0x00000006)
            @test layermask_signed(Float32, 8)[1] === UInt32(0x0000000E)
            @test layermask_signed(Float32, 2^31)[1] === UInt32(0xFFFFFFFE)

            # Float64 → UInt64
            @test layermask_signed(Float64, 1)[1] === UInt64(0x0000000000000000)
            @test layermask_signed(Float64, 2)[1] === UInt64(0x0000000000000002)
            @test layermask_signed(Float64, 4)[1] === UInt64(0x0000000000000006)
            @test layermask_signed(Float64, 8)[1] === UInt64(0x000000000000000E)
            @test layermask_signed(Float64, UInt64(2)^63)[1] === UInt64(0xFFFFFFFFFFFFFFFE)
        end
        @testset "second return value indicates overlap" begin
            # Float16 → UInt16
            @test layermask_signed(Float16, 1)[2] === false
            @test layermask_signed(Float16, 2)[2] === false
            @test layermask_signed(Float16, 4)[2] === false
            @test layermask_signed(Float16, 8)[2] === false
            @test layermask_signed(Float16, 2^5)[2] === false
            @test layermask_signed(Float16, 2^6)[2] === true
            @test layermask_signed(Float16, 2^15)[2] === true

            # Float32 → UInt32
            @test layermask_signed(Float32, 1)[2] === false
            @test layermask_signed(Float32, 2)[2] === false
            @test layermask_signed(Float32, 4)[2] === false
            @test layermask_signed(Float32, 8)[2] === false
            @test layermask_signed(Float32, 2^8)[2] === false
            @test layermask_signed(Float32, 2^9)[2] === true
            @test layermask_signed(Float32, 2^31)[2] === true

            # Float64 → UInt64
            @test layermask_signed(Float64, 1)[2] === false
            @test layermask_signed(Float64, 2)[2] === false
            @test layermask_signed(Float64, 4)[2] === false
            @test layermask_signed(Float64, 8)[2] === false
            @test layermask_signed(Float64, 2^11)[2] === false
            @test layermask_signed(Float64, 2^12)[2] === true
            @test layermask_signed(Float64, big(2)^63)[2] === true
        end
    end
end

nothing