# ZigguratTools
[![CI](https://github.com/npbarnes/ZigguratTools/actions/workflows/CI.yml/badge.svg)](https://github.com/npbarnes/ZigguratTools/actions/workflows/CI.yml)

This package will be a collection of tools for generating fast ziggurat-type random samplers for a large class of probability distributions. **Currently in development.**

## Goals

This package aims to make using the Ziggurat Method[^1] for random variate generation as simple as possible. Julia provides implementations of the ziggurat algorithm for normal and exponential distributions, but the algorithm could be adapted to a large class of distributions. The annoying part is generating the tables of values used by the algorithm. 

The table generation algorithm requires several inputs: pdf, inverse pdf, cdf, and mode. In addition, a fallback algorithm for the tail is needed (most likely using an inverse cdf). Having the user figure all that out is too much to ask, in my opinion. I want to automate it as much as possible. I plan to use root finding, autodifferentiation, and numerical integration to compute the inverse pdf, cdf, mode, and inverse cdf. Ideally, I'd like to be able to provide a pdf and get a sampler back that implements a ziggurat algorithm with performance similar to Julia's `randn` and `randexp` functions.  I expect to be able to achieve that high level of performance for sampling, but the ziggurat generation will be a potentially slow operation. Therefore, this algorithm will only make sense in contexts where many samples are needed from a fixed distribution.

At first, I am focusing on monotonic distributions with finite density. That includes functions that are not strictly monotonic. The ziggurat algorithm is usually applied to unimodal distributions by randomly selecting a sign, but I believe I can extend that to piecewise monotonic distributions using an alias table. In the future, I may also implement the Generalized Ziggurat Method from Jalavand and Charsooghi[^2] to support distributions with unbounded densities.

## Status
I would currently describe the status of this package with three words:
 * Incomplete
 * Buggy
 * Slow

But I am making progress on all three of those fronts.

I'm currently working on monotonic distributions. They will be the foundation of more complicated distributions, so it is important to get them right. Right now I can often make simple ziggurats and sample from them correctly, but there remain a lot of dangerous corner cases. For example, I started with a few manually written inverse pdfs, but I've found that floating point rounding can cause problems in specific circumstances. Counterintuitively, I think it's better to compute the generalized inverse using a bisection method since it can make certain guarantees that can't be made with floating point algebra. I would have liked to use Roots.jl or NonlinearSolve.jl for this problem, but they're inappropriate for this job because they return if they find an exact zero. In my case, having *an* x that satisfies the equation isn't enough. I need the *largest* x that satisfies the equation. That's a slightly different problem. Fortunately, a bisection search is not a complicated algorithm to implement.

## Installation
If it's not already obvious, this package isn't ready for widespread use, but you can play around with it if you'd like.

I intend to register v0.1 once I have most of the basic features working. Until then you can install it by tracking the main branch of this repo. Open a Julia REPL and type `]` to enter package mode. Then run
```julia
pkg> add https://github.com/npbarnes/ZigguratTools#main
```

## Examples
These examples appear to work right now. That doesn't mean that everything you try will work.

#### Distributions with Bounded Support
First define a pdf and its inverse. The pdf does not need to be normalized.
```julia-repl
julia> f(x) = 0 <= x <= 1 ? exp(-x^2) : zero(float(x))
julia> inv_f(x) = min(1.0, sqrt(-log(y)))
```
Then build the ziggurat
```julia-repl
julia> using ZigguratTools
julia> z = monotonic_ziggurat(10, 0.0, 1.0, f, inv_f)
```
`z` can now be used like any other sampler.
```julia-repl
julia> rand(z)
0.2977446038532221

julia> rand(z, 3)
3-element Vector{Float64}:
 0.07855690949819316      
 0.6499899606875755       
 0.7515162419547353       

julia> using Random

julia> a = Vector{Float64}(undef, 3)
...

julia> rand!(z, a)
3-element Vector{Float64}:
 0.18796323126287684
 0.21786647926588895
 0.0636838420645357

julia> rng = MersenneTwister(1234)
...

julia> rand(rng, z, 3)
3-element Vector{Float64}:
 0.6271541821806008
 0.6265042348161904
 0.005040015479099047
```

You can visualize the ziggurat layers
```julia-repl
julia> plotziggurat(z)
```
<img src="/assets/BoundedMonotonicZiggurat.svg" width=350/>
(right now Plots.jl is a dependency of ZigguratTools.jl just to define this function. Whenever I get around to it, I will make a plotting recipe with RecipesBase.jl. That should help with loading time.) 

You can also qualitatively validate the distribution with a histogram. Make sure everything is correctly normalized!
```julia-repl
julia> using Plots
julia> histogram(rand(z, 10^6), norm=:pdf)
julia> plot!(x -> f(x)/0.746824, color=:black, lw=3)
```
<img src="/assets/BoundedMonotonicZiggurat_Histogram.svg" width=350/>

[^1]: Marsaglia, G., & Tsang, W. W. (2000). The Ziggurat Method for Generating Random Variables. Journal of Statistical Software, 5(8), 1–7. https://doi.org/10.18637/jss.v005.i08
[^2]: Jalalvand, M., & Charsooghi, M. A. (2018). Generalized ziggurat algorithm for unimodal and unbounded probability density functions with Zest. arXiv preprint arXiv:1810.04744.
