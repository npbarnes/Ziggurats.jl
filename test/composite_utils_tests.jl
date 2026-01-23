@testset "CompositeZiggurat construction Utilities" begin
    @testset "get_subdomains" begin
        using Ziggurats: get_subdomains
        function no_discontinuity(x)
            2 - x
        end

        function left_closed(x)
            x <= 1 ? 0.0 : 1.0
        end

        function right_closed(x)
            x < 1 ? 0.0 : 1.0
        end

        function both_open(x)
            if x < 1
                return 0.0
            elseif x > 1
                return 1.0
            else
                return 0.5
            end
        end

        function nan_discontinuity(x)
            if x < 1
                return 0.0
            elseif x > 1
                return 1.0
            else
                return NaN
            end
        end

        funcs = [no_discontinuity, left_closed, right_closed, both_open, nan_discontinuity]

        @testset "Returns an n×2 Matrix for n subdomains" begin
            @testset for f in funcs
                @testset for n in 1:10
                    sd = get_subdomains(f, range(0, 2, length = (n+1)))
                    @test sd isa Matrix
                    @test size(sd) == (n, 2)
                end
            end
        end

        @testset "The first entry is the minimum of the domain" begin
            @testset for f in funcs
                d = regularize((0, 1, 2))
                m = get_subdomains(f, d)
                @test m[1, 1] == d[1]
            end
        end

        @testset "The last entry is the maximum of the domain" begin
            @testset for f in funcs
                d = regularize((0, 1, 2))
                m = get_subdomains(f, d)
                @test m[end, end] == d[end]
            end
        end

        @testset "Subdomains are split around a discontinuity so that it is not present in either half" begin
            @testset for f in funcs
                sd = get_subdomains(f, (0, 1, 2))
                @test f(sd[2, 1]) == 0.0
                @test f(sd[1, 2]) == 1.0
            end
        end

        @testset "When there is no discontinuity, the subdomains share the boundary point" begin
            sd = get_subdomains(no_discontinuity, (0, 1, 2))
            @test sd[2, 1] == sd[1, 2]
        end
    end
end

nothing
