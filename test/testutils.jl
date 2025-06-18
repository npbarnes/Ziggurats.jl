# These utilities are copied from:
# https://github.com/JuliaStats/Distributions.jl/blob/47c040beef8c61bad3e1eefa4fc8194e3a62b55a/test/testutils.jl
#
# I removed the functions/tests that I didn't need or want, and formatted them
# to match my codebase.

# Utilities to support the testing of distributions and samplers

using Distributions
using Random
using Printf: @printf
using Test: @test

# to workaround issues of Base.linspace
function _linspace(a::Float64, b::Float64, n::Int)
    intv = (b - a) / (n - 1)
    r = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        r[i] = a + (i - 1) * intv
    end
    r[n] = b
    return r
end

#################################################
#
#   Driver functions
#
#################################################

# testing the implementation of a continuous univariate distribution
#
function test_distr(
    distr::ContinuousUnivariateDistribution,
    n::Int;
    testquan::Bool = true,
    rng::AbstractRNG = MersenneTwister(123)
)
    test_range(distr)
    vs = get_evalsamples(distr, 0.01, 2000)

    test_support(distr, vs)
    test_evaluation(distr, vs, testquan)
    test_nonfinite(distr)

    if isa(distr, StudentizedRange)
        n = 2000 # must use fewer values due to performance
    end
    xs = test_samples(distr, n)
    xs = test_samples(distr, n; rng = rng)
    test_params(distr)
end

#################################################
#
#   Core testing functions
#
#################################################

#### Testing sampleable objects (samplers)

# for continuous samplers
#
function test_samples(
    s,    # the sampleable instance
    distr::ContinuousUnivariateDistribution,  # corresponding distribution
    n::Int = 10000;                                   # number of samples to generate
    nbins::Int = 50,                            # divide the main interval into nbins
    q::Float64 = 1.0e-6,                        # confidence interval, 1 - q as confidence
    verbose::Bool = false,                      # show intermediate info (for debugging)
    rng::Union{AbstractRNG,Missing} = missing
) # add an rng?

    # The basic idea
    # ------------------
    #   Generate n samples, divide the interval [q00, q99] into nbins bins, where
    #   q01 and q99 are 0.01 and 0.99 quantiles, and count the numbers of samples
    #   falling into each bin. For each bin, we will compute a confidence interval
    #   of the number, and see whether the actual number is in this interval.
    #
    #   If the distribution has a bounded range, it also checks whether
    #   the samples are all within this range.
    #
    #   By setting a small q, we ensure that failure of the tests rarely
    #   happen in practice.
    #

    if verbose
        println("test_samples on $(typeof(s))")
    end

    if !(n > 1)
        error("The number of samples must be greater than 1.")
    end
    if !(nbins > 1)
        error("The number of bins must be greater than 1.")
    end
    if !(0.0 < q < 0.1)
        error("The value of q must be within the open interval (0.0, 0.1).")
    end

    # determine the range of values to examine
    vmin = minimum(distr)
    vmax = maximum(distr)

    local rmin::Float64
    local rmax::Float64
    if applicable(quantile, distr, 0.5)
        rmin = quantile(distr, 0.01)
        rmax = quantile(distr, 0.99)
    elseif isfinite(vmin) && isfinite(vmax)
        rmin = vmin
        rmax = vmax
    end
    edges = _linspace(rmin, rmax, nbins + 1)

    # determine confidence intervals for counts:
    # with probability q, the count will be out of this interval.
    #
    clb = Vector{Int}(undef, nbins)
    cub = Vector{Int}(undef, nbins)
    cdfs = map(Base.Fix1(cdf, distr), edges)

    for i in 1:nbins
        pi = cdfs[i + 1] - cdfs[i]
        bp = Binomial(n, pi)
        clb[i] = floor(Int, quantile(bp, q / 2))
        cub[i] = ceil(Int, cquantile(bp, q / 2))
        @assert cub[i] >= clb[i]
    end

    # generate samples using RNG passed or default RNG
    # we also check reproducibility
    if rng === missing
        Random.seed!(1234)
        samples = rand(s, n)
        Random.seed!(1234)
        samples2 = rand(s, n)
    else
        rng2 = deepcopy(rng)
        samples = rand(rng, s, n)
        samples2 = rand(rng2, s, n)
    end
    @test length(samples) == n
    @test samples2 == samples

    if isa(distr, StudentizedRange)
        samples[isnan.(samples)] .= 0.0 # Underlying implementation in Rmath can't handle very low values.
    end

    # check whether all samples are in the valid range
    for i in 1:n
        @inbounds si = samples[i]
        if !(vmin <= si <= vmax)
            error("Sample value out of valid range.")
        end
    end

    # get counts
    cnts = fit(Histogram, samples, edges; closed = :right).weights
    @assert length(cnts) == nbins

    # check the counts
    for i in 1:nbins
        if verbose
            @printf(
                "[%.4f, %.4f) ==> (%d, %d): %d\n",
                edges[i],
                edges[i + 1],
                clb[i],
                cub[i],
                cnts[i]
            )
        end
        if !(clb[i] <= cnts[i] <= cub[i])
            error("The counts are out of the confidence interval.")
        end
    end

    #@test eltype(samples) === eltype(distr)
    @test !any(isnan, samples)
    @test !any(isinf, samples) skip = eltype(s)==Float16

    return samples
end

function test_samples(
    distr::ContinuousUnivariateDistribution,
    n::Int;
    nbins::Int = 50,
    q::Float64 = 1.0e-6,
    verbose::Bool = false,
    rng = missing
)
    test_samples(distr, distr, n; nbins = nbins, q = q, verbose = verbose, rng = rng)
end

#### Testing range & support methods

function test_range(d::UnivariateDistribution)
    vmin = minimum(d)
    vmax = maximum(d)
    @test vmin <= vmax

    is_lb = islowerbounded(d)
    is_ub = isupperbounded(d)

    @test isfinite(vmin) == is_lb
    @test isfinite(vmax) == is_ub
    @test isbounded(d) == (is_lb && is_ub)
end

function get_evalsamples(d::ContinuousUnivariateDistribution, q::Float64, n::Int)
    # samples for testing evaluation functions (even spacing)

    lv = quantile(d, q / 2)
    hv = cquantile(d, q / 2)
    @assert isfinite(lv) && isfinite(hv) && lv <= hv
    return _linspace(lv, hv, n)
end

function test_support(d::UnivariateDistribution, vs::AbstractVector)
    for v in vs
        @test insupport(d, v)
    end
    @test all(insupport(d, vs))

    if islowerbounded(d)
        @test isfinite(minimum(d))
        @test insupport(d, minimum(d))
        @test !insupport(d, minimum(d) - 1)
    end
    if isupperbounded(d)
        @test isfinite(maximum(d))
        @test insupport(d, maximum(d))
        @test !insupport(d, maximum(d) + 1)
    end

    @test isbounded(d) == (isupperbounded(d) && islowerbounded(d))

    # Test the `Base.in` or `∈` operator
    # The `support` function is buggy for unbounded `DiscreteUnivariateDistribution`s
    if isbounded(d) || isa(d, ContinuousUnivariateDistribution)
        s = support(d)
        for v in vs
            @test v ∈ s
        end

        if islowerbounded(d)
            @test minimum(d) ∈ s
            @test (minimum(d) - 1) ∉ s
        end
        if isupperbounded(d)
            @test maximum(d) ∈ s
            @test (maximum(d) + 1) ∉ s
        end
    end
end

#### Testing evaluation methods

function test_evaluation(
    d::ContinuousUnivariateDistribution,
    vs::AbstractVector,
    testquan::Bool = true
)
    nv = length(vs)
    p = Vector{Float64}(undef, nv)
    c = Vector{Float64}(undef, nv)
    cc = Vector{Float64}(undef, nv)
    lp = Vector{Float64}(undef, nv)
    lc = Vector{Float64}(undef, nv)
    lcc = Vector{Float64}(undef, nv)

    for (i, v) in enumerate(vs)
        if !isa(d, StudentizedRange)
            p[i] = pdf(d, v)
            lp[i] = logpdf(d, v)
            @assert p[i] >= 0.0
        end

        c[i] = cdf(d, v)
        cc[i] = ccdf(d, v)
        lc[i] = logcdf(d, v)
        lcc[i] = logccdf(d, v)

        @assert (i == 1 || c[i] >= c[i - 1])

        @test isapprox(c[i] + cc[i], 1.0, atol = 1.0e-12)
        if !isa(d, StudentizedRange)
            @test isapprox(lp[i], log(p[i]), atol = 1.0e-12)
        end
        @test isapprox(lc[i], log(c[i]), atol = 1.0e-12)
        @test isapprox(lcc[i], log(cc[i]), atol = 1.0e-12)

        if testquan
            # TODO: remove this line when we have more accurate implementation
            # of quantile for InverseGaussian and StudentizedRange
            qtol = isa(d, InverseGaussian) ? 1.0e-4 : 1.0e-10
            qtol = isa(d, StudentizedRange) ? 1.0e-5 : qtol
            if p[i] > 1.0e-6
                @test isapprox(quantile(d, c[i]), v, atol = qtol * (abs(v) + 1.0))
                @test isapprox(cquantile(d, cc[i]), v, atol = qtol * (abs(v) + 1.0))
                @test isapprox(invlogcdf(d, lc[i]), v, atol = qtol * (abs(v) + 1.0))
                @test isapprox(invlogccdf(d, lcc[i]), v, atol = qtol * (abs(v) + 1.0))
            end
        end
    end

    if !isa(d, StudentizedRange)
        # check: pdf should be the derivative of cdf
        for i in 2:(nv - 1)
            if p[i] > 1.0e-6
                v = vs[i]
                ap = (cdf(d, v + 1.0e-6) - cdf(d, v - 1.0e-6)) / (2.0e-6)
                @test isapprox(p[i], ap, atol = p[i] * 1.0e-3)
            end
        end
    end

    # consistency of scalar-based and vectorized evaluation
    if !isa(d, StudentizedRange)
        @test Base.Fix1(pdf, d).(vs) ≈ p
        @test Base.Fix1(logpdf, d).(vs) ≈ lp
    end

    @test Base.Fix1(cdf, d).(vs) ≈ c
    @test Base.Fix1(ccdf, d).(vs) ≈ cc

    @test Base.Fix1(logcdf, d).(vs) ≈ lc
    @test Base.Fix1(logccdf, d).(vs) ≈ lcc
end

function test_nonfinite(distr::UnivariateDistribution)
    # non-finite bounds
    @test iszero(@inferred(cdf(distr, -Inf)))
    @test isone(@inferred(cdf(distr, Inf)))
    @test isone(@inferred(ccdf(distr, -Inf)))
    @test iszero(@inferred(ccdf(distr, Inf)))
    @test @inferred(logcdf(distr, -Inf)) == -Inf
    @test iszero(@inferred(logcdf(distr, Inf)))
    @test iszero(@inferred(logccdf(distr, -Inf)))
    @test @inferred(logccdf(distr, Inf)) == -Inf

    # NaN
    for f in (cdf, ccdf, logcdf, logccdf)
        if distr isa NoncentralT
            # broken in StatsFuns/R
            @test_broken isnan(f(distr, NaN))
        else
            @test isnan(f(distr, NaN))
        end
    end
end

#### Testing statistics methods

function test_params(d::Distribution)
    # simply test that params returns something sufficient to
    # reconstruct d
    D = typeof(d)
    pars = params(d)
    d_new = D(pars...)
    @test d_new == d
    @test d_new == deepcopy(d)
end

function test_params(d::Truncated)
    # simply test that params returns something sufficient to
    # reconstruct d
    d_unt = d.untruncated
    D = typeof(d_unt)
    pars = params(d_unt)
    d_new = truncated(D(pars...), d.lower, d.upper)
    @test d_new == d
    @test d == deepcopy(d)
end
