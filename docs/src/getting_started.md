# Getting Started

## Installation
You can install the package with the Julia package manager Pkg:
```
# Press ']' to enter the Pkg REPL mode.
pkg> add PolyaGammaHybridSamplers
```
or:
```julia
julia> using Pkg
julia> Pkg.add("PolyaGammaHybridSamplers")
```

## Sampling

Start by including the package:
```julia
julia> using PolyaGammaHybridSamplers
```
Then create a sampler object:
```julia
julia> s = PolyaGammaHybridSampler(5, 4.0)
```
The `rand` function can be used to draw samples from the sampler object `s`:
```julia
julia> rand(s, 3)
```
The default sampling method is `HYBRID`, which dymanically selects between `DEVROYE`, `SADDLEPOINT`, and `NORMALAPPROX` depending on the parameters (see the README for more details). A sampler can be forced to always use a certain method by including it as an argument, e.g.:
```julia
julia> s = PolyaGammaHybridSampler(5, 4.0, SADDLEPOINT)
```
