using ZigguratTools, Test

include("./testutils.jl")

dist = SteppedExponential()

# This distribution, with this N causes an error if we use a naive search algorithm
@test_nowarn UnboundedDecreasingZiggurat(
    x->pdf(dist,x),
    y->ipdf_right(dist,y),
    x->ccdf(dist,x),
    mode(dist),
    256,
    x->sampler(truncated(dist, lower=x))
)
