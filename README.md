# ZigguratTools
[![CI](https://github.com/npbarnes/ZigguratTools/actions/workflows/CI.yml/badge.svg)](https://github.com/npbarnes/ZigguratTools/actions/workflows/CI.yml)

This package contains routines for automatically constructing high-performance non-uniform
random number generators for a wide variety of univariate probability distributions using a
variation on the Marsaglia & Tsang Ziggurat Algorithm[^1].

The `ziggurat` function provides the primary interface. In most cases you can pass in the
`pdf` and a list of points that divides the domain into monotonic segments.
```julia-repl
julia> using ZigguratTools, Plots
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
In some cases, the monotonic segments can be found automatically using autodifferentiation
and root finding. The public non-exported function `monotonic_segments` does exactly that.
```julia-repl
julia> import ZigguratTools: monotonic_segments
julia> f(x) = exp(-x^2/2) * (1+sin(3x)^2) * (1+cos(5x)^2);
julia> d = monotonic_segments(f, (-3,3));
julia> z = ziggurat(f, d);
julia> histogram(rand(z, 10^6); normalize=:pdf, fillrange=0);
julia> plot!(x->f(x)/5.62607; linewidth=2, color=:black)
```
<img src="/assets/multimodal.svg" width=350/><br/>
Discontinuous functions are also okay. They do not need to be divided into subdomains at the
discontinuity as long as it's monotonic.
```julia-repl
julia> z = ziggurat(x -> x>=1 ? 2-x : 10-x, (0,2));
julia> histogram(rand(z, 10^5), normalize=:pdf);
julia> plot!(x -> (x>=1 ? 2-x : 10-x)/10; linewidth=2, color=:black)
```
<img src="/assets/discontinuity.svg" width=350/><br/>

# How it Works
For monotonic distributions, the algorithm is essentially the same as Marsaglia & Tsang 2000[^1]
with a few minor technical differences. The most important difference from the classic
algorithm is that the inverse pdf function, tail area function, and fallback algorithm don't
need to be explicitly provided by the programmer. This makes the ziggurat method much more
ergonomic. By default, the `ipdf` is computed using a root finding algorithm (via `Roots.jl`),
and the `tailarea` function is computed using Gauss-Kronrod quadrature (via `QuadGK.jl`).
The fallback algorithm is inverse transform sampling over the renormalized `tailarea` function,
where the inverse is computed similarly to the `ipdf`.

For piecewise monotonic functions the same process is used to make a ziggurat for each
monotonic segment. Then the relative area under the curve of each segment is computed by
quadrature. To sample a point, a ziggurat is selected using an Alias Table of the relative
areas (`AliasTables.jl`) and then that ziggurat is sampled just as in the monotonic case. The
classic Marsaglia & Tsang algorithm is highly optimized for the symmetric unimodal case, but
those optimizations are not yet available in ZigguratTools.jl. This algorithm is more flexible
at the expense of some performance.

# Performance Tips
Generally, you can expect high performance, but there are a couple of things you should be
aware of to get the best possible performance. 

First of all, the most optimized samplers are available when the number of layers is a power
of two with `N <= 4096` for `Float64`, `N <= 512` for `Float32`, and `N <= 64` for `Float16`.
This applies to ziggurats on both bounded and unbounded domains.

Secondly, on unbounded domains, the default fallback algorithm is usually by far the slowest
part. Conventional wisdom with the ziggurat algorithm is that the performance of the fallback
algorithm is not particularly important because it is called very rarely. However, in this
case, the fallback can be exceptionally slow. It often *does* have a meaningful impact on the
average time to generate a sample. Fortunately, there are several ways to mitigate the problem:

1. **Use more layers**. The more layers there are, the more rarely the fallback is required,
the smaller its impact on overall performance. This is probably a good choice in most
situations unless your application is memory constrained or sensitive to the maximum possible
latency.

2. **Use a bounded domain**. Ziggurats on bounded domains don't use a fallback at all, so
you can just pick an outer limit that's beyond the majority of the probability mass. This
could cut off a portion of the tail, but you can probably always find a cutoff point that is
far enough down the tail that it makes no material difference. For example, with the
exponential distribution, `exp(-1000) == 0` in floating point, so chances are using a domain
of `(0,1000)` would incur no loss of precision.

3. **Provide a function for the tail area**. The default fallback is so slow because it's a
numerical inverse of a numerical integral. Replacing the tail area function with a direct
implementation will greatly improve the performance of the fallback as well.

4. **Override the fallback directly**. If you happen to know a modestly performant way to
sample points in the tail, you can use it. You're on your own for an algorithm, though.

The next section demonstrates each of these approaches.

# Benchmarks
As an initial example here is the exponential distribution
```julia-repl
julia> using BenchmarkTools
julia> z = ziggurat(x->exp(-x), (0, Inf));
julia> @benchmark rand($z)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations per sample.
 Range (min … max):   3.930 ns … 693.475 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):      4.660 ns               ┊ GC (median):    0.00%
 Time  (mean ± σ):   68.164 ns ±  94.573 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

 █                  ▂▇▁                 ▁▄                    ▁
 █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▇███▇▆▄▁▃▃▃▁▁▃▁▁▁▄▃▅▆███▇▅▄▃▁▁▁▃▁▁▁▁▁▁▁▁▅▆ █
 3.93 ns       Histogram: log(frequency) by time       415 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
If you are familiar with `BenchmarkTools.jl` you may have been conditioned to focus your
attention on the minimum time (at least, I was). In that case, `ZigguratTools` compares well
against `Random.jl`'s hand optimized `randexp`
```julia-repl
julia> using Random
julia> @benchmark randexp()
BenchmarkTools.Trial: 10000 samples with 1000 evaluations per sample.
 Range (min … max):  4.851 ns … 13.852 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     5.311 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   5.341 ns ±  0.333 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

            ▁▃▁▅▂█▅▆█▃▆▁▄▂                                    
  ▂▂▂▃▃▃▄▄▆▆██████████████▇█▆▆▅▄▄▃▃▃▂▂▂▂▂▂▂▂▂▁▁▂▂▂▁▂▂▂▁▁▂▂▂▂ ▄
  4.85 ns        Histogram: frequency by time        6.35 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
However, the mean and maximum times reveal a discrepancy. In this case, the minimum and
median times are representative of the fast path which executes the vast majority of the time,
but the maximum and mean times reveal that the slower fallback method for generating samples
in the tail is much slower than the specialized method `randexp()`. This can be mitigated
using the methods described in the previous section.
```julia-repl
julia> z = ziggurat(x->exp(-x), (0, Inf), 4096);
julia> @benchmark rand($z)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations per sample.
 Range (min … max):  3.739 ns … 276.213 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     3.829 ns               ┊ GC (median):    0.00%
 Time  (mean ± σ):   6.708 ns ±  19.660 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

 █                                                         ▁ ▁
 █▄▄▄▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▄▄█ █
 3.74 ns      Histogram: log(frequency) by time       140 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> z = ziggurat(x->exp(-x), (0, 1000));
julia> @benchmark rand($z)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations per sample.
 Range (min … max):  4.100 ns … 23.897 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     4.639 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   4.658 ns ±  0.355 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

            ▁▃▁▅▂█▅▆█▃▆▁▄▂                                    
  ▂▂▂▃▃▃▄▄▆▆██████████████▇█▆▆▅▄▄▃▃▃▂▂▂▂▂▂▂▂▂▁▁▂▂▂▁▂▂▂▁▁▂▂▂▂ ▄
  4.85 ns        Histogram: frequency by time        6.35 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> z = ziggurat(x->exp(-x), (0, Inf); tailarea=x->exp(-x));
julia> @benchmark rand($z)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations per sample.
 Range (min … max):  4.049 ns … 37.065 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     4.640 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   5.524 ns ±  1.553 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▄█▇▁                                                      
  ▃▆████▄▃▂▂▂▂▂▂▁▂▂▂▂▃▅██▆▄▃▂▂▂▂▂▁▂▁▂▁▂▂▂▃▃▃▂▂▂▂▂▂▂▁▁▁▁▂▁▂▂▂ ▃
  4.05 ns        Histogram: frequency by time        11.1 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> z = ziggurat(x->exp(-x), (0, Inf); fallback = (rng, a) -> a - log1p(-rand(rng)));
julia> @benchmark rand($z)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations per sample.
 Range (min … max):  3.970 ns … 22.527 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     4.509 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   4.526 ns ±  0.356 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

                    ▁ ▅▃▆▆█▅▃▇▇▇█▄▃▆▅▅▃▁                      
  ▁▁▁▁▁▁▁▁▂▂▃▃▃▃▅▆▆██▅██████████████████▅█▇▆▆▅▃▃▃▃▃▂▂▁▂▂▁▁▁▁ ▄
  3.97 ns        Histogram: frequency by time        5.04 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
Now we have attained a similar performance as `randexp()`. This is to be expected since
`randexp` uses a very similar ziggurat method internally.

Piecewise-monotonic distributions are implemented as a collection of monotonic ziggurats.
Each of the unbounded segments can be optimized in any of the same ways. However, the gains
are not as substantial since selecting a ziggurat adds constant overhead.
```julia-repl
julia> z = ziggurat(x -> exp(-x^2/2), (-40, 0, 40));
julia> @benchmark rand($z)
BenchmarkTools.Trial: 10000 samples with 997 evaluations per sample.
 Range (min … max):  18.987 ns … 106.098 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     19.819 ns               ┊ GC (median):    0.00%
 Time  (mean ± σ):   19.921 ns ±   1.267 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

         ▁▄███▇▆▁                                               
  ▂▂▂▂▃▄▇████████▇▆▄▃▃▃▃▂▂▂▂▂▂▂▂▂▁▁▂▁▁▂▂▁▂▂▁▂▂▁▁▁▁▁▁▁▁▁▁▂▁▂▂▂▂ ▃
  19 ns           Histogram: frequency by time         23.4 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
Generally speaking, these results are pretty fast, but it's not as impressive compared to a
specialized algorithm like Julia's `randn()`.

That's because the classic Marsaglia & Tsang algorithm is highly optimized for the symmetric
unimodal case. `ZigguratTools` doesn't currently implement those optimizations. However, in
exchange, we have much greater flexibility to sample from a wide variety of distributions.
Adding more segments comes with no additional overhead. Here's a benchmark of one of the
examples above with 18 segments:
```julia-repl
julia> f(x) = exp(-x^2/2) * (1+sin(3x)^2) * (1+cos(5x)^2);
julia> d = monotonic_segments(f, (-3,3));
julia> z = ziggurat(f, d);
julia> @benchmark rand($z)
BenchmarkTools.Trial: 10000 samples with 998 evaluations per sample.
 Range (min … max):  18.357 ns … 55.000 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     19.499 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   19.694 ns ±  1.499 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

      ▃▇█▆▂                                                    
  ▂▂▃▆█████▆▄▃▂▂▂▂▂▁▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂ ▃
  18.4 ns         Histogram: frequency by time        28.7 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
It has essentially identical performance as the normal distribution example.

# Contributing
 - **Star this repo!** This is a hobby project. If you star the repo it lets me know that
 people are interested and that I should keep working on it.

 - **Raise an issue.** Any input helps. If you think something could be better, let me know!
 If you have a feature request, let me know! If you find a bug then please post it in the
 issue tracker!

 - **Send a pull request.** If you can contribute code or documentation (even small edits)
 then please do so. I want this project to be as good as it can be, so I'm happy to take
 outside contributions.

[^1]: Marsaglia, G., & Tsang, W. W. (2000). The Ziggurat Method for Generating Random Variables.
Journal of Statistical Software, 5(8), 1–7. https://doi.org/10.18637/jss.v005.i08
