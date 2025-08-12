@testset "Test Test Distributions" begin
    # Test the test-distributions
    @testset "Stepped Exponential Distributions" begin
        test_distr(SteppedExponential(), 10^6)
        test_distr(SteppedExponential(0.5), 10^6)
        test_distr(SteppedExponential(1.5), 10^6)
    end

    @testset "Doorstop Distribution" begin
        test_distr(Doorstop(-1, 1, 3), 10^6)
        test_distr(Doorstop(-1, -1, 3), 10^6)
        test_distr(Doorstop(-1, 1, 1), 10^6)
    end
end
