@testset "layer_bits()" begin
    import Ziggurats: layer_bits, layer_bits_signed, layermask, layermask_signed

    int64gen = Data.Integers{UInt64}()
    int32gen = Data.Integers{UInt32}()
    int16gen = Data.Integers{UInt16}()

    @testset "When the Layer Mask only uses lower bits, there is no overlap and no rearrangement" begin
        @testset "Float64" begin
            @testset "unsigned" begin
                @testset for nbits in 0:12
                    nlayers = UInt64(2)^nbits
                    LM = layermask(Float64, nlayers)
                    @check function nooverlap(r = int64gen)
                        layer_bits(Float64, LM, r) == r & LM
                    end
                end
            end

            @testset "signed" begin
                @testset for nbits in 0:11
                    nlayers = UInt64(2)^nbits
                    LM = layermask_signed(Float64, nlayers)
                    @check function nooverlap_signed(r = int64gen)
                        layer_bits_signed(Float64, LM, r) == (r & LM) >>> 1
                    end
                end
            end
        end

        @testset "Float32" begin
            @testset "unsigned" begin
                @testset for nbits in 0:9
                    nlayers = UInt32(2)^nbits
                    LM = layermask(Float32, nlayers)
                    @check function nooverlap(r = int32gen)
                        layer_bits(Float32, LM, r) == r & LM
                    end
                end
            end

            @testset "signed" begin
                @testset for nbits in 0:8
                    nlayers = UInt32(2)^nbits
                    LM = layermask_signed(Float32, nlayers)
                    @check function nooverlap_signed(r = int32gen)
                        layer_bits_signed(Float32, LM, r) == (r & LM) >>> 1
                    end
                end
            end
        end

        @testset "Float16" begin
            @testset "unsigned" begin
                @testset for nbits in 0:6
                    nlayers = UInt16(2)^nbits
                    LM = layermask(Float16, nlayers)
                    @check function nooverlap(r = int16gen)
                        layer_bits(Float16, LM, r) == r & LM
                    end
                end
            end

            @testset "signed" begin
                @testset for nbits in 0:5
                    nlayers = UInt16(2)^nbits
                    LM = layermask_signed(Float16, nlayers)
                    @check function nooverlap_signed(r = int16gen)
                        layer_bits_signed(Float16, LM, r) == (r & LM) >>> 1
                    end
                end
            end
        end
    end

    @testset "The non-overlapped bits are preserved as the upper bits of the result" begin
        @testset "Float64" begin
            @testset "unsigned" begin
                @testset for nbits in 13:64
                    nlayers = big(2)^nbits # 2^64 needs BigInt to be represented
                    LM = layermask(Float64, nlayers)
                    @check function preserve_nonoverlap(r = int64gen)
                        nonoverlap_mask = UInt64(2)^12 - UInt64(1)
                        nonoverlapbits = r & nonoverlap_mask
                        amountofoverlap = nbits - 12

                        re = layer_bits(Float64, LM, r)

                        (re >>> amountofoverlap) == nonoverlapbits
                    end
                end
            end

            @testset "signed" begin
                @testset for nbits in 12:63
                    nlayers = UInt64(2)^nbits
                    LM = layermask_signed(Float64, nlayers)
                    @check function preserve_nonoverlap_signed(r = int64gen)
                        nonoverlap_mask = (UInt64(2)^11 - UInt64(1)) << 1
                        nonoverlapbits = (r & nonoverlap_mask) >>> 1
                        amountofoverlap = nbits - 11

                        re = layer_bits_signed(Float64, LM, r)

                        (re >>> amountofoverlap) == nonoverlapbits
                    end
                end
            end
        end

        @testset "Float32" begin
            @testset "unsigned" begin
                @testset for nbits in 10:32
                    nlayers = UInt64(2)^nbits # 2^32 needs UInt64 to be represented
                    LM = layermask(Float32, nlayers)
                    @check function preserve_nonoverlap(r = int32gen)
                        nonoverlap_mask = UInt32(2)^9 - UInt32(1)
                        nonoverlapbits = r & nonoverlap_mask
                        amountofoverlap = nbits - 9

                        re = layer_bits(Float32, LM, r)

                        (re >>> amountofoverlap) == nonoverlapbits
                    end
                end
            end

            @testset "signed" begin
                @testset for nbits in 9:31
                    nlayers = UInt32(2)^nbits
                    LM = layermask_signed(Float32, nlayers)
                    @check function preserve_nonoverlap_signed(r = int32gen)
                        nonoverlap_mask = (UInt32(2)^8 - UInt32(1)) << 1
                        nonoverlapbits = (r & nonoverlap_mask) >>> 1
                        amountofoverlap = nbits - 8

                        re = layer_bits_signed(Float32, LM, r)

                        (re >>> amountofoverlap) == nonoverlapbits
                    end
                end
            end
        end

        @testset "Float16" begin
            @testset "unsigned" begin
                @testset for nbits in 7:16
                    nlayers = UInt32(2)^nbits # 2^16 needs UInt32 to be represented
                    LM = layermask(Float16, nlayers)
                    @check function preserve_nonoverlap(r = int16gen)
                        nonoverlap_mask = (UInt16(2)^6 - UInt16(1))
                        nonoverlapbits = r & nonoverlap_mask
                        amountofoverlap = nbits - 6

                        re = layer_bits(Float16, LM, r)

                        (re >>> amountofoverlap) == nonoverlapbits
                    end
                end
            end

            @testset "signed" begin
                @testset for nbits in 6:15
                    nlayers = UInt16(2)^nbits
                    LM = layermask_signed(Float16, nlayers)
                    @check function preserve_nonoverlap_signed(r = int16gen)
                        nonoverlap_mask = (UInt16(2)^5 - UInt16(1)) << 1
                        nonoverlapbits = (r & nonoverlap_mask) >>> 1
                        amountofoverlap = nbits - 5

                        re = layer_bits_signed(Float16, LM, r)

                        (re >>> amountofoverlap) == nonoverlapbits
                    end
                end
            end
        end
    end

    @testset "The overlapped bits are preserved as the lower bits of the result" begin
        @testset "Float64" begin
            @testset "unsigned" begin
                @testset for nbits in 13:64
                    nlayers = big(2)^nbits # 2^64 needs BigInt to be represented
                    LM = layermask(Float64, nlayers)
                    amountofoverlap = nbits - 12
                    lowbit_mask = UInt64(2^amountofoverlap - 1)
                    @check function perserve_overlap(r = int64gen)
                        re = layer_bits(Float64, LM, r)
                        re & lowbit_mask == (r >>> 12) & lowbit_mask
                    end
                end
            end

            @testset "signed" begin
                @testset for nbits in 12:63
                    nlayers = UInt64(2)^nbits
                    LM = layermask_signed(Float64, nlayers)
                    amountofoverlap = nbits - 11
                    lowbit_mask = UInt64(2^amountofoverlap - 1)
                    @check function perserve_overlap_signed(r = int64gen)
                        re = layer_bits_signed(Float64, LM, r)
                        re & lowbit_mask == (r >>> 12) & lowbit_mask
                    end
                end
            end
        end

        @testset "Float32" begin
            @testset "unsigned" begin
                @testset for nbits in 10:32
                    nlayers = UInt64(2)^nbits # 2^32 needs UInt64 to be represented
                    LM = layermask(Float32, nlayers)
                    amountofoverlap = nbits - 9
                    lowbit_mask = UInt32(2^amountofoverlap - 1)
                    @check function perserve_overlap(r = int32gen)
                        re = layer_bits(Float32, LM, r)
                        re & lowbit_mask == (r >>> 9) & lowbit_mask
                    end
                end
            end
            @testset "signed" begin
                @testset for nbits in 9:31
                    nlayers = UInt32(2)^nbits
                    LM = layermask_signed(Float32, nlayers)
                    amountofoverlap = nbits - 8
                    lowbit_mask = UInt32(2^amountofoverlap - 1)
                    @check function perserve_overlap_signed(r = int32gen)
                        re = layer_bits_signed(Float32, LM, r)
                        re & lowbit_mask == (r >>> 9) & lowbit_mask
                    end
                end
            end
        end

        @testset "Float16" begin
            @testset "unsigned" begin
                @testset for nbits in 7:16
                    nlayers = UInt32(2)^nbits # 2^16 needs UInt32 to be represented
                    LM = layermask(Float16, nlayers)
                    amountofoverlap = nbits - 6
                    lowbit_mask = UInt16(2^amountofoverlap - 1)
                    @check function perserve_overlap(r = int16gen)
                        re = layer_bits(Float16, LM, r)
                        re & lowbit_mask == (r >>> 6) & lowbit_mask
                    end
                end
            end
            @testset "signed" begin
                @testset for nbits in 6:15
                    nlayers = UInt16(2)^nbits
                    LM = layermask_signed(Float16, nlayers)
                    amountofoverlap = nbits - 5
                    lowbit_mask = UInt16(2^amountofoverlap - 1)
                    @check function perserve_overlap_signed(r = int16gen)
                        re = layer_bits_signed(Float16, LM, r)
                        re & lowbit_mask == (r >>> 6) & lowbit_mask
                    end
                end
            end
        end
    end

    @testset "for 2^m layers, the m+1 and higher bits are unset" begin
        @testset "Float64" begin
            @testset "unsinged" begin
                @testset for nbits in 13:64
                    nlayers = big(2)^nbits # 2^64 needs BigInt to be represented
                    LM = layermask(Float64, nlayers)
                    @check function highbits_unset(r = int64gen)
                        re = layer_bits(Float64, LM, r)
                        (re & ~UInt64(nlayers - 1)) == 0
                    end
                end
            end
            @testset "signed" begin
                @testset for nbits in 12:63
                    nlayers = UInt64(2)^nbits
                    LM = layermask_signed(Float64, nlayers)
                    @check function highbits_unset_signed(r = int64gen)
                        re = layer_bits_signed(Float64, LM, r)
                        (re & ~(nlayers - 1)) == 0
                    end
                end
            end
        end
        @testset "Float32" begin
            @testset "unsigned" begin
                @testset for nbits in 10:32
                    nlayers = UInt64(2)^nbits # 2^32 needs UInt64 to be represented
                    LM = layermask(Float32, nlayers)
                    @check function highbits_unset(r = int32gen)
                        re = layer_bits(Float32, LM, r)
                        (re & ~UInt32(nlayers - 1)) == 0
                    end
                end
            end

            @testset "signed" begin
                @testset for nbits in 9:31
                    nlayers = UInt32(2)^nbits
                    LM = layermask_signed(Float32, nlayers)
                    @check function highbits_unset_signed(r = int32gen)
                        re = layer_bits_signed(Float32, LM, r)
                        (re & ~(nlayers - 1)) == 0
                    end
                end
            end
        end
        @testset "Float16" begin
            @testset "unsigned" begin
                @testset for nbits in 7:16
                    nlayers = UInt32(2)^nbits # 2^16 needs UInt32 to be represented
                    LM = layermask(Float16, nlayers)
                    @check function highbits_unset(r = int16gen)
                        re = layer_bits(Float16, LM, r)
                        (re & ~UInt16(nlayers - 1)) == 0
                    end
                end
            end
            @testset "signed" begin
                @testset for nbits in 6:15
                    nlayers = UInt16(2)^nbits
                    LM = layermask_signed(Float16, nlayers)
                    @check function highbits_unset_signed(r = int16gen)
                        re = layer_bits_signed(Float16, LM, r)
                        (re & ~(nlayers - 1)) == 0
                    end
                end
            end
        end
    end
end

nothing
