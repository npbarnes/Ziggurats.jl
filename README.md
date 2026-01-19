# Ziggurats.jl
[![CI](https://github.com/npbarnes/Ziggurats.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/npbarnes/Ziggurats.jl/actions/workflows/CI.yml)

This package contains routines for automatically constructing high-performance non-uniform
random number generators for a wide variety of user defined univariate probability distributions using a
variation on the Marsaglia & Tsang Ziggurat Algorithm[^1].

It's designed as a black-box algorithm, in most cases all the programmer needs to provide is
a probability density function (pdf) and a list of extremal points of the pdf that divide the
domain into subdomains where the pdf is monotonic.

The Ziggurat Algorithm is designed for very high marginal sampling performance at the expense
of a longer setup time. It's appopriate for situations where you need a large number of samples
from a fixed distribution with a known pdf.

# Installation
Ziggurats.jl is a registered Julia package that can be installed using the Julia package manager.
E.g., in the REPL type
```julia-repl
julia> using Pkg; Pkg.add("Ziggurats")
```

# Basic Usage
The `ziggurat` function provides the primary interface. In most cases you can pass in the
`pdf` and a list of extremal points. That is, the boundaries of the domain and each local
minimum and maximum.
```julia-repl
julia> using Ziggurats, Plots
julia> z = ziggurat(x -> exp(-x^2/2), (-Inf, 0, Inf));
julia> histogram(rand(z, 10^5); normalize=:pdf);
julia> plot!(x->1/√(2π) * exp(-x^2/2); linewidth=2, color=:black)
```
```julia-repl
julia> z = ziggurat(x -> abs(cos(x)), (-π, -π/2, 0, π/2, π));
julia> histogram(rand(z, 10^6); normalize=:pdf)
julia> plot!(x -> 1/4 * abs(cos(x)); linewidth=2, color=:black)
```
<img src="/assets/normal.svg" width=350/> <img src="/assets/cos.svg" width=350/> <br/>

Discontinuous functions are okay. Points of discontinuity do not need to be specified unless
it's also an extremal point.
```julia-repl
julia> z = ziggurat(x -> x>=1 ? 2-x : 10-x, (0,2));
julia> histogram(rand(z, 10^5), normalize=:pdf);
julia> plot!(x -> (x>=1 ? 2-x : 10-x)/10; linewidth=2, color=:black)
```
<img src="/assets/discontinuity.svg" width=350/><br/>

For symmetric and unimodal pdfs (bell-shaped), a more optimized sampler can be
generated using the `BellZiggurat` constructor. Pass in the pdf function and a pair of points
that are the mode and one boundary of the domain. Either boundary is okay.
```julia-repl
julia> z = BellZiggurat(x -> exp(-x^2/2), (0, Inf));
# Or
julia> z = BellZiggurat(x -> exp(-x^2/2), (-Inf, 0));
```

The available constructors for ziggurats are: `ziggurat`, `BellZiggurat`, `monotonic_ziggurat`,
`BoundedZiggurat`, `UnboundedZiggurat`, and `CompositeZiggurat`. Each of them has a detailed
docstring where you can find information about usage and limitations.

# How it Works
For monotonic and bell-shaped distributions, the algorithm is an easy to use, robust, high performance,
implementation of the Marsaglia & Tsang 2000[^1] Ziggurat Algorithm. Unlike most implementations,
Ziggurats.jl is not confined to a small number of distributions. It's is a 'black-box' algorithm
that automates the ziggurat construction producing fast samplers for a wide variety of
user-specified distributions. Usually only a pdf function and a list of extremal points of the pdf
is required. Other libraries with similar
functionality require the programmer to provide an inverse pdf function, a cumulative distribution
function, etc.. Ziggurats.jl by default generates those auxilurary functions using robust
numerical methods. This makes the ziggurat method much more ergonomic to use.

By default, the `ipdf` is computed
using a root finding algorithm (via `Roots.jl`), and the `tailarea` function is computed using
Gauss-Kronrod quadrature (via `QuadGK.jl`). The fallback algorithm for samples in the tail of
unbounded distributions is inverse transform sampling over the renormalized `tailarea` function,
where the inverse is computed similarly to the `ipdf`.

Ziggurats.jl also extends the ziggurat algorithm to distributions with peicewise monotonic pdfs.
The extension is relatively simple, ziggurats are generated for each monotonic subdomain
(which are separated by extremal points), and
their areas are computed by quadtrature. To take a sample, first select a subdomain with the appropriate probability
(proportional to area, via `AliasTables.jl`) and then sample within that subdomain using the
corresponding ziggurat.

# Performance Tips
Generally, you can expect high performance, but there are a few things you should be
aware of to get the best possible performance.

First of all, if your pdf is known to be symmetric and unimodal (bell-shaped) use `BellZiggurat`
instead of `ziggurat`. The specialized algorithm can be several times faster than the generic
algorithm.

Second, bulk generation is faster than generating samples one at a time, so if you're filling
an array with random numbers use `rand(z, dims)`, or `rand!(a, z)` instead of `[rand(z) for _ in 1:10^6]`,
or
```julia
for i in eachindex(a)
    a[i] = rand(z)
end
```

Third, using more layers in the ziggurat may be faster because it reduces the rejection rate.
Any number of layers is possible,
but the fastest sampling algorithms are available when the number of layers is a power of two.
I recommend benchmarking your application with different powers of two numbers of layers up to a maximum of 2^12 = 4096
for Float64 domains, 2^9 = 512 for Float32 domains and 2^6 = 64
for Float16 domains. For techincal reasons, the quality of the statistical distribution can degrade if you
use a power of two larger than those limits. The default number of layers is 2^8 = 256 or 2^6=64 for Float16.
That's generally a good trade off between speed and memory size.

Finally, on unbounded domains, the default fallback algorithm is usually by far the slowest
part. Conventional wisdom with the ziggurat algorithm is that the performance of the fallback
algorithm is not particularly important because it is called very rarely. However, in this
case, the fallback can be exceptionally slow. It often *does* have a meaningful impact on the
average time to generate a sample. Fortunately, there are several ways to mitigate the problem:

1. **Use more layers**. The more layers there are, the more rarely the fallback is required,
the smaller its impact on overall performance. This is probably a good choice in most
situations unless your application is memory constrained (including cache memory) or sensitive to the maximum possible
latency. I recommend benchmarking your application with different powers of
two numbers of layers up to the limits mentioned above.

2. **Use a bounded domain**. Ziggurats on bounded domains don't use a fallback at all, so
you can just pick an outer limit that's beyond the majority of the probability mass. This
could cut off a portion of the tail, but you can probably always find a cutoff point that is
far enough down the tail that it makes no material difference. For example, with the
exponential distribution, `exp(-1000) == 0` in floating point, so chances are using a domain
of `(0,1000)` would incur no loss of precision.

3. **Provide a function for the tail area**. The default fallback is so slow because it's a
numerical inverse of a numerical integral. Replacing the tail area function with a direct
implementation by passing a `tailarea` keyword argument will greatly improve the performance of the fallback as well.

4. **Override the fallback directly**. If you happen to know a modestly performant way to
sample points in the tail, you can use it by passing a `fallback` keyword argument. See the
constructor docstrings for details on usage.

# Contributing
 - **Star this repo!** This is a hobby project. If you star the repo it lets me know that
 people are interested and that I should keep working on it.

 - **Start a discussion.** If you need help, or have a question, click on the discussion tab
 on the GitHub repo and make a post. I'll do my best to help out.

 - **Raise an issue.** Any input helps. If you think something could be better, let me know!
 If you have a feature request, let me know! If you find a bug then please post it in the
 issue tracker!

 - **Send a pull request.** If you can contribute code or documentation (even small edits)
 then please do so. I want this project to be as good as it can be, so I'm happy to take
 outside contributions.

[^1]: Marsaglia, G., & Tsang, W. W. (2000). The Ziggurat Method for Generating Random Variables.
Journal of Statistical Software, 5(8), 1–7. https://doi.org/10.18637/jss.v005.i08
