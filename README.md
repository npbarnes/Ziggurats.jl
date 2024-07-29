# ZigguratTools
[![Runtests](https://github.com/npbarnes/ZigguratTools/actions/workflows/Runtests.yml/badge.svg)](https://github.com/npbarnes/ZigguratTools/actions/workflows/Runtests.yml)

This package aims to be a collection of tools for generating fast ziggurat-type random samplers for arbitrary distributions. **Currently in development.**

## Goals

This package aims to make using the Ziggurat Method[^1] for random variate generation as simple as possible. Julia provides implementations of the ziggurat algorithm for normal and exponential distributions, but the algorithm could be applied to a large class of distributions. The annoying part is generating the tables of values used by the algorithm. The table generation algorithm requires several inputs: pdf, inverse pdf, cdf, and mode. In addition, a fallback algorithm for the tail is needed (most likely using an inverse cdf). Having the user figure all that out is too much to ask, in my opinion. I want to automate it as much as possible. Ideally, I'd like to be able to provide a pdf and get a sampler back that implements a ziggurat algorithm with performance similar to Julia's `randn` and `randexp` functions. The plan is to use root finding, autodifferentiation, and numerical integration to compute the inverse pdf, cdf, mode, and inverse cdf. There can be edge cases where, for example, QuadGK doesn't work very well, but in general, I expect this approach to work for most pdf's. 

At first, I will focus on monotonic and unimodal distributions with finite density. In the future, I may also implement the Generalized Ziggurat Method of Jalavand and Charsooghi[^2] to support distributions with unbounded densities.

[^1]: Marsaglia, G., & Tsang, W. W. (2000). The Ziggurat Method for Generating Random Variables. Journal of Statistical Software, 5(8), 1–7. https://doi.org/10.18637/jss.v005.i08
[^2]: Jalalvand, M., & Charsooghi, M. A. (2018). Generalized ziggurat algorithm for unimodal and unbounded probability density functions with Zest. arXiv preprint arXiv:1810.04744.
