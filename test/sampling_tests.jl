using ZigguratTools, Test
using Distributions
using Random
using StatsBase

function test_samples(s::AbstractZiggurat,    # the sampleable instance
    distr::ContinuousUnivariateDistribution,  # corresponding distribution
    n::Int;                                   # number of samples to generate
    nbins::Int=50,                            # divide the main interval into nbins
    q::Float64=1.0e-6,                        # confidence interval, 1 - q as confidence
    verbose::Bool=false,                      # show intermediate info (for debugging)
    rng::Union{AbstractRNG,Missing}=missing) # add an rng?

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

    verbose && println("test_samples on $(typeof(s))")

    n > 1 || error("The number of samples must be greater than 1.")
    nbins > 1 || error("The number of bins must be greater than 1.")
    0.0 < q < 0.1 || error("The value of q must be within the open interval (0.0, 0.1).")

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
    edges = range(rmin, rmax, nbins + 1)

    # determine confidence intervals for counts:
    # with probability q, the count will be out of this interval.
    #
    clb = Vector{Int}(undef, nbins)
    cub = Vector{Int}(undef, nbins)
    cdfs = map(Base.Fix1(cdf, distr), edges)

    for i = 1:nbins
        pi = cdfs[i+1] - cdfs[i]
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
    for i = 1:n
        @inbounds si = samples[i]
        vmin <= si <= vmax ||
            error("Sample value out of valid range.")
    end

    # get counts
    cnts = fit(Histogram, samples, edges; closed=:right).weights
    @assert length(cnts) == nbins

    # check the counts
    for i = 1:nbins
        if verbose
            println("[$(edges[i]), $(edges[i+1])) ==> ($(clb[i]), $(cub[i])): $(cnts[i])")
        end
        clb[i] <= cnts[i] <= cub[i] ||
            error("The counts are out of the confidence interval.")
    end
    return samples
end

function testsampling(dist, z)
    values = test_samples(z, dist, 10000)

    @test eltype(values) === eltype(dist)
    @test !any(isnan, values)
    @test !any(isinf, values)
end

@testset "Normal (x>=0)" begin
    basedist = Normal()

    # Because of the choice of UnboundedDecreasingZiggurat, ipdf_right, and ccdf,
    # this ziggurat will actually be sampling from truncated(Normal(), lower=0.0).
    # A small number N is chosen so that the fallback branch gets chosen with high
    # probability and its confidence intervals can be tested.
    z = UnboundedDecreasingZiggurat(
        x->pdf(basedist,x),
        y->ipdf_right(basedist,y),
        x->ccdf(basedist,x),
        mode(basedist),
        2,
        x->sampler(truncated(basedist, lower=x))
    )

    testsampling(truncated(basedist, lower=mode(basedist)), z)
end

@testset "Normal (x<=0)" begin
    basedist = Normal()

    # Because of the choice of a AbstractUnboundedMonotonicZiggurat, ipdf_left, and cdf,
    # this ziggurat will actually be sampling from truncated(Normal(), lower=0.0).
    # A small number N is chosen so that the fallback branch gets chosen with high
    # probability and its confidence intervals can be tested.
    z = UnboundedIncreasingZiggurat(
        x->pdf(basedist,x),
        y->ipdf_left(basedist,y),
        x->cdf(basedist,x),
        mode(basedist),
        2,
        x->sampler(truncated(basedist, upper=x))
    )

    testsampling(truncated(basedist, upper=mode(basedist)), z)
end

@testset "Exponential" begin
    basedist = Exponential()

    z = UnboundedDecreasingZiggurat(
        x->pdf(basedist,x),
        y->ipdf_right(basedist,y),
        x->ccdf(basedist,x),
        mode(basedist),
        2,
        x->sampler(truncated(basedist, lower=x))
    )

    testsampling(basedist, z)
end

@testset "SteppedExponential" begin
    include("./testutils.jl")
    
    basedist = SteppedExponential()

    z = UnboundedDecreasingZiggurat(
        x->pdf(basedist,x),
        y->ipdf_right(basedist,y),
        x->ccdf(basedist,x),
        mode(basedist),
        2,
        x->sampler(truncated(basedist, lower=x))
    )

    testsampling(basedist, z)
end