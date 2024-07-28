using ZigguratTools, Test

include("./testutils.jl")

dist = SteppedExponential()

# This distribution, with this N causes an error if we use a naive search
# algorithm. For an continuous pdf, `zig.y = pdf.(zig.x)`, but when the pdf is
# discontinuous, that need not be the case. We only require `zig.x = ipdf.(zig.y)`, for a
# generalized inverse. Therefore, there are multiple valid y's for each x at the
# discontinuity. The ziggurat search algorithm needs to be aware of this. With
# SteppedExponential and 256 layers, the correct first layer needs to have
# pdf(prevfloat(8.0)) < y < pdf(8.0)
@test_nowarn monotonic_ziggurat(
    256,
    mode(dist),
    x->ccdf(dist,x),
    x->pdf(dist,x),
    y->ipdf_right(dist,y),
    x->sampler(truncated(dist, lower=x))
)
