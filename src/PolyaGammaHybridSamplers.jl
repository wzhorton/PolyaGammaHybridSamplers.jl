module PolyaGammaHybridSamplers

using Distributions
using Random
using SpecialFunctions

export PolyaGammaHybridSampler
export PGSamplingMethod, HYBRID, DEVROYE, SADDLEPOINT, NORMALAPPROX
export rand_pggammasum

@enum PGSamplingMethod HYBRID DEVROYE SADDLEPOINT NORMALAPPROX

include("rand_pgdevroye.jl")
include("rand_pgsaddlepoint.jl")
include("rand_pgnormalapprox.jl")
include("rand_pggammasum.jl")
include("pghybrid.jl")

end
