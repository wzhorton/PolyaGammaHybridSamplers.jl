module PolyaGammaHybridSamplers

using Distributions
using Random

export PolyaGammaHybridSampler, PGSamplingMethod

include("rand_pgdevroye.jl")
include("rand_pgsaddlepoint.jl")
include("rand_pgnormalapprox.jl")
include("pghybrid.jl")

end
