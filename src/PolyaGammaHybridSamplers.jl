module PolyaGammaHybridSamplers

using Distributions
using Random

export PolyaGammaHybridSampler
export PGSamplingMethod, HYBRID, DEVROYE, SADDLEPOINT, NORMALAPPROX

include("rand_pgdevroye.jl")
include("rand_pgsaddlepoint.jl")
include("rand_pgnormalapprox.jl")
include("pghybrid.jl")

end
