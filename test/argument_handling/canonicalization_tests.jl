using Test

import Ziggurats:
    CanonicalArray,
    RollingPairs,
    ManifestlyListOfListOfScalar,
    ManifestlyListOfScalar,
    ListOfScalar,
    ListOfListOfScalar,
    Scalar,
    SomethingElse
import Ziggurats: canonicalize, best_float_type, _best_float_type, classify, _classify, islistofscalar, _islistofscalar

@testset "RollingPairs Iterator Tests" begin
    @testset "Edge cases: 0 and 1 elements" begin
        # Zero elements
        r0 = RollingPairs([])
        @test collect(r0) == []
        @test length(r0) == 0

        # One element must throw an explicit descriptive error
        err_msg = "a RollingPairs iterator cannot be formed from an iterator with just one element. There must be zero, 2, or more elements."
        @test_throws ErrorException(err_msg) RollingPairs([1])
    end

    @testset "Standard operation: 2 and multiple elements" begin
        # Two elements (1 pair)
        r2 = RollingPairs([10, 20])
        @test collect(r2) == [(10, 20)]
        @test length(r2) == 1
        @test eltype(r2) == Tuple{Int,Int}

        # Multiple elements
        r3 = RollingPairs([1, 2, 3, 4])
        @test collect(r3) == [(1, 2), (2, 3), (3, 4)]
        @test length(r3) == 3
        @test eltype(r3) == Tuple{Int,Int}

        # Mixed types
        r4 = RollingPairs((1, 2, 3.0))
        @test eltype(r4) == Tuple{Real,Real}
        r5 = RollingPairs(Real[1, 2, 3.0])
        @test eltype(r5) == Tuple{Real,Real}
    end

    @testset "Custom HasLength() iterator support" begin
        gen = (x * 2 for x in 1:4)
        r_gen = RollingPairs(gen)
        @test collect(r_gen) == [(2, 4), (4, 6), (6, 8)]
        @test length(r_gen) == 3
    end
end

@testset "Classification and Type Helpers" begin
    @testset "use type info to shortcut classification" begin
        # ManifestlyListOfScalar
        MLS = [[1, 2, 3], Number[1, 2, 3], (1, 2, 3), (), Float64[]]

        for mls in MLS
            @test mls isa ManifestlyListOfScalar
            @test @inferred(classify(mls)) isa ListOfScalar
            @test _classify(mls) isa ListOfScalar
            @test islistofscalar(mls)
            @test _islistofscalar(mls)
        end

        # ManifestlyListOfListOfScalar
        MLLS = [
            [[1, 2], [3, 4]],
            [(1, 2), (3, 4)],
            ([1, 2], [3, 4]),
            ((1, 2), (3, 4)),
            [Number[1, 2], Number[3, 4]],
            [(1, 2), (3, 4, 5)],
            Union{AbstractVector{<:Number},Tuple{Number,Number}}[[1, 2], (3, 4)],
            ([1, 2], (3, 4)),
            [[1, 2, 3]],
            ([1, 2, 3],),
            [(1, 2, 3)],
            ((1, 2, 3),),
            [Float64[]],
            [()],
            (Float64[],),
            ((),)
        ]

        @testset for mlls in MLLS
            @test mlls isa ManifestlyListOfListOfScalar
            @test @inferred(classify(mlls)) isa ListOfListOfScalar
            @test _classify(mlls) isa ListOfListOfScalar
            @test !islistofscalar(mlls)
            @test !_islistofscalar(mlls)
            @test all(islistofscalar, mlls)
            @test all(_islistofscalar, mlls)
        end
    end

    @testset "classify / islistofscalar" begin
        # Scalars
        scalars = (42, 3.14)
        @testset for s in scalars
            @test classify(s) isa Scalar
            @test classify(s) == _classify(s)
        end

        # List of Scalars
        listofscalar = ([1, 2, 3], Any[1, 2, 3], (1.0, 2.0, 3.0), (i^2 for i in Any[1, 2, 3]))
        @testset for ls in listofscalar
            @test classify(ls) isa ListOfScalar
            @test classify(ls) == _classify(ls)
        end

        # List of List of Scalars
        listoflistofscalar = (
            [[1, 2, 3]],
            [(1, 2, 3)],
            ([1, 2, 3],),
            ((1, 2, 3),),
            [[1, 2], [3, 4]],
            Any[Any[1, 2], Any[3, 4]],
            ((1, 2), (3, 4)),
            ([1, 2], (3, 4)),
            ((i+j for j in Any[0, 1]) for i in Any[1, 2, 3])
        )
        @testset for lls in listoflistofscalar
            @test classify(lls) isa ListOfListOfScalar
            @test classify(lls) == _classify(lls)
        end

        # Something Else (Strings, unsupported types, wrong nesting)
        somethingelse = (
            "not a list",
            :notiterable,
            [:a, :b],
            [Ref(:a), Ref(:b)],
            (1, (2, 3)),
            (((1, 2), (3, 4)), ((5, 6),)),
            Any[1, Any[2, 3]],
            Any[Any[Any[1, 2], Any[3, 4]], Any[Any[5, 6]]]
        )
        @testset for se in somethingelse
            @test classify(se) isa SomethingElse
            @test classify(se) == _classify(se)
        end
    end

    @testset "best_float_type promotion" begin
        inferrable_cases = [
            (1, 2, 3) => Float64,
            [1, 2, 3] => Float64,
            Float32[1, 2, 3] => Float32,
            (1, 2, 3.0f0) => Float32,
            ((1, 2), (3, 4)) => Float64,
            [[1, 2], [3, 4]] => Float64,
            ((1, 2), (3, 4.0f0)) => Float32,
            ([1, 2], Float32[3, 4]) => Float32,
            ([1, 2], [3, 4]) => Float64,
            (Int[1, 2], Float32[3, 4]) => Float32,
            ([1, 2], (3, 4.0f0)) => Float32,
            (1.0f0, Float16(2.0)) => Float32,
            (1, π) => Float64,
            (1.0f0, π) => Float32
        ]

        noninferrable_cases = [
            Any[1, 2, 3] => Float64,
            Any[1, 2, 3.0f0] => Float32,
            (Number[1, 2], [3, 4]) => Float64,
            Vector{Real}[[1, 2], Float32[3, 4]] => Float32,
            [Any[1, 2], Any[3, 4.0f0]] => Float32,
            Real[1.0f0, Float16(2.0)] => Float32,
            Real[1, π] => Float64,
            Real[1.0f0, π] => Float32
        ]

        @testset for tc in inferrable_cases
            d = tc.first
            T = tc.second
            @test @inferred(best_float_type(d)) == T
            @test best_float_type(d) == _best_float_type(d)

            # canonicalize is type inferred whenever best_float_type is type inferred.
            @test @inferred(canonicalize(d)) isa CanonicalArray{T}
        end

        @testset for tc in noninferrable_cases
            d = tc.first
            T = tc.second
            @test_throws ErrorException @inferred(best_float_type(d))
            @test best_float_type(d) == T
            @test best_float_type(d) == _best_float_type(d)
        end
    end
end

@testset "Validation and Error Messages" begin
    @testset "ListOfScalar Validation Failures" begin
        # Contiguous domain too short (< 2 elements)
        @test_throws ErrorException("a contiguous domain needs at least two distinct points to mark the boundaries. Got [1].") canonicalize([1])
        @test_throws ErrorException("a contiguous domain needs at least two distinct points to mark the boundaries. Got [3].") canonicalize([
            [1, 2],
            [3]
        ])

        # Contiguous domain not sorted
        @test_throws ErrorException("the critical points that make up a contiguous domain must be in increasing order. Got [2, 1].") canonicalize([
            2,
            1
        ])
        @test_throws ErrorException("the critical points that make up a contiguous domain must be in increasing order. Got [2, 1].") canonicalize([
            [2, 1],
            [2, 3]
        ])

        # Duplicates found
        @test_throws ErrorException("duplicate values found in a contiguous domain = [1, 1, 2].") canonicalize([
            1,
            1,
            2
        ])

        # Invalid element types inside purported scalar list
        @test_throws ErrorException("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got [\"a\", \"b\"].") canonicalize([
            "a",
            "b"
        ])
    end

    @testset "ListOfListOfScalar Validation Failures (Overlaps)" begin
        # Overlapping subdomains
        overlapping1 = [[1, 3], [2, 4]]
        err_msg1 = "the contiguous domains that make up `domain` must be non-overlapping. Got $overlapping1."
        @test_throws ErrorException(err_msg1) canonicalize(overlapping1)

        overlapping2 = [[2, 4], [1, 3]]
        err_msg2 = "the contiguous domains that make up `domain` must be non-overlapping. Got $overlapping2."
        @test_throws ErrorException(err_msg2) canonicalize(overlapping2)
    end

    @testset "Top-level Domain Validation Failures" begin
        # Passing a scalar directly
        @test_throws ErrorException("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got 42.") canonicalize(42)

        # Passing an unsupported type (SomethingElse)
        @test_throws ErrorException("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got [:a, :b].") canonicalize([
            :a,
            :b
        ])
    end
end

@testset "Canonicalization with Supported Types" begin
    @testset "Flat numbers (Vectors and Tuples)" begin
        @test canonicalize([1, 2, 3]) == [(1.0, 2.0), (2.0, 3.0)]
        @test canonicalize((1, 2, 3, 4)) == [(1.0, 2.0), (2.0, 3.0), (3.0, 4.0)]
    end

    @testset "Views" begin
        base_vec = [10, 20, 30, 40, 50]
        sub_view = @view base_vec[2:4] # [20, 30, 40]
        @test canonicalize(sub_view) == [(20.0, 30.0), (30.0, 40.0)]
    end

    @testset "Iterators with length()" begin
        @test canonicalize((i+j for j in Any[0, 1]) for i in Any[1, 2, 3]) == [(1, 2), (2, 3), (3, 4)]
    end

    @testset "Nested collections (Vectors of Vectors / Tuples)" begin
        # Vector of vectors
        @test canonicalize([[1, 2], [3, 4, 5]]) == [(1.0, 2.0), (3.0, 4.0), (4.0, 5.0)]
        @test canonicalize([[1, 2], (3, 4, 5)]) == [(1.0, 2.0), (3.0, 4.0), (4.0, 5.0)]
        @test canonicalize(((1, 2, 3), [4, 5])) == [(1.0, 2.0), (2.0, 3.0), (4.0, 5.0)]
        # Each contiguous domain must in sorted order, but several contiguous domains need
        # not be in any particular order
        @test canonicalize([[3, 4, 5], [1, 2]]) == [(1.0, 2.0), (3.0, 4.0), (4.0, 5.0)]
        @test canonicalize([[1, 2], (3, 4, 5), [2, 3]]) == [(1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.0)]
        @test canonicalize(([4, 5], (1, 2, 3))) == [(1.0, 2.0), (2.0, 3.0), (4.0, 5.0)]
    end
end

@testset "Wrong Nesting and Descriptive Errors" begin
    @testset "Mixed scalar and list at top level" begin
        bad_mix = [(1, 2), 3]
        @test_throws ErrorException("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got $bad_mix.") canonicalize(bad_mix)
    end

    @testset "Over-nested lists (lists of lists of lists)" begin
        over_nested = [[[1, 2], [3, 4]], [5, 6]]
        @test_throws ErrorException("a domain must be either a iterator of Numbers or an iterator of iterators of Numbers. Got $over_nested.") canonicalize(over_nested)
    end
end

@testset "canonicalization is idempotent" begin
    # calling canonicalize() on a CanonicalArray will return it imediately
    # calling canonicalize() on an array that is in the canonical form will be wrapped and returned.
    a = canonicalize([1, 2, 3, 4])
    @test canonicalize(a) === a
    @test canonicalize(a.arr).arr === a.arr
    @test @allocated(canonicalize(a)) == 0
    @test @allocated(canonicalize(a.arr)) == 0
end

@testset "canonicalization is type-inferred for tuples and some concrete array types" begin
    inferrable = [
        (1, 2, 3),
        [1, 2, 3],
        (1.0, 2.0),
        (1.0f0, 2.0f0),
        (1.0f0, 2),
        Float32[1, 2, 3],
        [[1, 2], [3, 4]],
        [Float32[1, 2], Float32[3, 4]],
        [(1, 2), (3, 4)],
        ([1, 2], [3, 4]),
        (Float32[1, 2], Float16[3, 4]),
        ((1, 2), (3, 4)),
        [[1, 2], (3, 4)],
        ([1, 2], (3, 4)),
        [[1, 2, 3]],
        ([1, 2, 3],),
        [(1, 2, 3)],
        ((1, 2, 3),)
    ]
end
