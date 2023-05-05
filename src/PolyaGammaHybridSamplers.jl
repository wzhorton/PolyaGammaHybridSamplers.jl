module PolyaGammaHybridSamplers

using Distributions

export PolyaGammaHybridSampler

include("rand_pgdevroye.jl")
include("rand_pgsaddlepoint.jl")
include("rand_pgnormalapprox.jl")
include("pghybrid.jl")

end
