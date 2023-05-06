# `PolyaGammaHybridSamplers.jl`

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wzhorton.github.io/PolyaGammaHybridSamplers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wzhorton.github.io/PolyaGammaHybridSamplers.jl/dev/)
[![Build Status](https://github.com/wzhorton/PolyaGammaHybridSamplers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wzhorton/PolyaGammaHybridSamplers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/wzhorton/PolyaGammaHybridSamplers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/wzhorton/PolyaGammaHybridSamplers.jl)


`PolyaGammaHybridSamplers.jl` is a Julia package that implements a hybrid sampling approach for the Pólya-Gamma distribution. The Pólya-Gamma distribution is a continuous distribution on the positive real number line and is most often used for latent variable augmentation in models like logistic regression or logit stick-breaking processes. It is parameterized by two values, $b > 0$ and $z\in\mathbb{R}$, and is often assigned the distributional notation $PG(b,z)$. An excellent tutorial on using the Pólya-Gamma distribution can be found at this [blog post](https://gregorygundersen.com/blog/2019/09/20/polya-gamma/) by Gregory Gunderson.

# Background Details

The class of Pólya-Gamma distributions was first presented in a paper by [Polson et al. (2013)](http://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001). There no closed-form expression for the $PG(b,z)$ density function $f$, but is rather defined through its Laplace transform: 

$$
\mathcal{L}_f(t) = \frac{\cosh^b(\frac{z}{2})}{\cosh^b\left(\sqrt{\frac{z^2/2 + t}{2}}\right)}
$$

This can be derived from an alternative definition where the $PG(b,z)$ distribution is represented as an infinite mixture of gamma variables. The authors comment that a naive approach to sampling could be to generate atruncated number of gamma values and combine them using the mixture formula, but they warn that this approximation is both inefficient and potentially dangerous. 

In order to develop an effecient yet exact sampling method, the authors drew connections to a subclass of $J^*$ distributions discussed in [Devroye (2009)](https://www.sciencedirect.com/science/article/pii/S0167715209002867). The ultimate conclusion is an extremely efficient rejection sampler, called the *Devroye method*, that is used to simulate from a $PG(1,z)$ distribution. If $b$ is an integer, then the $PG(b,z)$ distribution can be drawn by summing $b$ independent draws from the $PG(1,z)$ distribution.

The following year, a technical report written by [Windle et al. (2014)](https://arxiv.org/abs/1405.0506) was published on arXiv detailing an alternative sampling technique and an approximate technique. The alternative technique requires $b\in[1,4]$, but allows for otherwise non-integer $b$ values. The approximate algorithm, called the *saddle-point sampler*, involves approximating the density using the moment generating function and then constructing a smart rejection sampler. The saddle-point sampler is both fast and increasingly accurate for larger values of $b$. The report concludes with a performance comparison between the *naive gamma mixture*, the *Devroye method*, the *alternative technique*, the *saddle-point sampler*, and a normal approximation. They subsequently suggest the following hybrid sampler scheme to balance accuracy and efficiency:

| Sampler | Ideal Circumstances |
| --- | --- |
| Devroye Method | $b = 1,2$ |
| Alternative Technique | $3 \le b \le 13$ and non-integer $1 < b < 3$ |
| Saddle-point Sampler| $13 < b < 170$ |
| Normal Approximation | $b \ge 170$ |

# Previous Code Implementations

The authors of [Polson et al. (2013)](http://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001) published the [BayesLogit](https://cran.r-project.org/web/packages/BayesLogit/index.html) R package with C++ code implementing the Devroye sampling method, later expanding it to include the techniques presented in [Windle et al. (2014)](https://arxiv.org/abs/1405.0506). The source code is available at the [BayesLogit GitHub site](https://github.com/jwindle/BayesLogit). However, it is important to note that, as of version 2.1, inline comments for the main sampling function, `rpg_hybrid`, state that the *alternative technique* code needs further review and/or development ([source](https://github.com/jwindle/BayesLogit/blob/master/src/polyagamma_wrapper.cpp), lines 138-141) and instead the *naive gamma mixture* is used for small, non-integer $b$ values.

Interest in the Pólya-Gamma distribution within the Julia community appears first in a [GitHub issues thread](https://github.com/JuliaStats/Distributions.jl/issues/685) started by a user called `currymj` under the main [Distributions.jl](https://github.com/JuliaStats/Distributions.jl/) repository. Ultimately, `currymj` developed a separate package called [PolyaGammaDistribution.jl](https://github.com/currymj/PolyaGammaDistribution.jl) that was essentially a direct translation of the [BayesLogit](https://cran.r-project.org/web/packages/BayesLogit/index.html) *Devroye method* code for integer $b$ parameters.

However, Julia would eventually update to v1.0 and begin requiring a different package structure. When asked about future compatibility, a collaborator, user `maximerischard`, indicated in a [GitHub issues thread](https://github.com/currymj/PolyaGammaDistribution.jl/issues/3) that there were no plans to upgrade the package. The user, `igutierrezm`, who started the thread would then go on to develop a Julia v1.0 compatible package, [PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl). This new package implemented the *Devroye method* similar to its predecessor, but using the new [sampler framework](https://juliastats.org/Distributions.jl/stable/extends/#Create-New-Samplers-and-Distributions), which is a more limited implementation compared to the full distribution requirements of [Distributions.jl](https://github.com/JuliaStats/Distributions.jl/), but is also more fitting for the Pólya-Gamma distribution.

Over time, a handful of variations came about, including a [PolyaGammaDistribution.jl fork](https://github.com/currymj/PolyaGammaDistribution.jl/forks) by user `yunzli` ([link](https://github.com/yunzli/PolyaGammaDistribution.jl)) and a rewritten *Devroye method* sampler in [AugmentedGaussianProcesses.jl](https://github.com/theogf/AugmentedGaussianProcesses.jl) by user `theogf`. Perhaps the most recent discourse comes from another [GitHub issues thread](https://github.com/JuliaStats/Distributions.jl/issues/1440), where user `theogf` asked whether their *Devroye method* code should be merged into [Distributions.jl](https://github.com/JuliaStats/Distributions.jl/). The issue was ultimately never resolved due to conflicts between the MIT and GPL-3.0 licenses.

# Scope of `PolyaGammaHybridSamplers.jl`

Nearly all existing native Julia implemetations of the Pólya-Gamma distribution focus exclusively on the *Devroye method* sampler (except perhaps for this [PolyaGammaSamplers.jl branch](https://github.com/igutierrezm/PolyaGammaSamplers.jl/tree/sp-approximation)). It is an exact sampler and is quite efficient for small values of $b$, but results from [Windle et al. (2014)](https://arxiv.org/abs/1405.0506) suggest that it is potentially orders of magnitude slower than the *saddle-point sampler* approximation for larger values of $b$. Additionally, none of the packages considered a normal approximation for large $b$, which is extremely efficient. In short, there is no hybrid sampler available in Julia like there is in the original [BayesLogit](https://cran.r-project.org/web/packages/BayesLogit/index.html) R package. As the name implies, the `PolyaGammaHybridSamplers.jl` package aims to provide one. The hybrid sampler implemented here is similar to the one given by [Windle et al. (2014)](https://arxiv.org/abs/1405.0506): 

| Sampler | Condition |
| --- | --- |
| Devroye Method | $b = 1,2,...,13$ |
| Saddle-point Sampler| $b = 13,..., 170$ |
| Normal Approximation | $b \ge 170$ |

To generate a random variate, create a sampler object and then use `rand` (vectorized versions of `rand` and `rand!` are also available):

```julia
pg = PolyaGammaHybridSampler(4, 1.5)
data = rand(pg)
```

--**The current plan is for the hybrid sampler to require integers for the $b$ parameter**.-- However, the underlying *normal approximation* and *saddle-point sampler* routines do support non-integer inputs. Additionally, the type of sampler can controlled by including, for example, `SADDLEPOINT` as an argument: 

```julia
pg = PolyaGammaHybridSampler(4, 1.5, SADDLEPOINT)
data = rand(pg)
```
Other valid options include `HYBRID`, `DEVROYE`, and `NORMALAPPROX`. But be careful as no warning about the approximation quality or efficiency will be given.

Perhaps in the future the *approximate technique* or the *gamma mixture* sampler will be implemented, maybe even with some clever multiple dispatching. However, this seems unlikely given that the vast majority of applications using the Pólya-Gamma distribution involve integer $b$, not to mention that the efficiency benefits of the *alternative technique* are pretty small and the fact that truncating the gamma mixture is "dangerous" according to [Polson et al. (2013)](http://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001).


# License Issue and Merging with Distributions.jl

This package is effectively a translation of the original [Polson et al. (2013)](http://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001) C++ code in the  package [BayesLogit](https://cran.r-project.org/web/packages/BayesLogit/index.html), including the later developments described in [Windle et al. (2014)](https://arxiv.org/abs/1405.0506). Some structural inspiration is also taken from the Julia packages [PolyaGammaDistribution.jl](https://github.com/currymj/PolyaGammaDistribution.jl) (`currymj`), [PolyaGammaSamplers.jl](https://github.com/igutierrezm/PolyaGammaSamplers.jl) (`igutierrezm`), and [AugmentedGaussianProcesses.jl](https://github.com/theogf/AugmentedGaussianProcesses.jl) (`theogf`). All of these packages are licensed under the GPL-3.0 license, except for [AugmentedGaussianProcesses.jl](https://github.com/theogf/AugmentedGaussianProcesses.jl) which is licensed under the [MIT "Expat" license](https://github.com/theogf/AugmentedGaussianProcesses.jl/blob/master/LICENSE.md). Therefore, the `PolyaGammaHybridSamplers.jl` package is necessarily licensed under GPL-3.0. The question of whether this code can merged into [Distributions.jl](https://github.com/JuliaStats/Distributions.jl/) has likely been answered in this [GitHub issues thread](https://github.com/JuliaStats/Distributions.jl/issues/1440) started by `theogf`: the MIT license of Julia's base environment conflicts with the GPL-3.0 requirements, so probably not.
