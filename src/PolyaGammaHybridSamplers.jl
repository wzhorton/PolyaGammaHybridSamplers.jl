module PolyaGammaHybridSamplers

using Distributions
using Random
using SpecialFunctions
using StatsFuns

export PolyaGammaHybridSampler
export PGSamplingMethod, DEVROYE, SADDLEPOINT, NORMALAPPROX, GAMMASUM, DEVROYEPLUSGAMMASUM 
export rand_pggammasum, rand_pgdevroye, rand_pgsaddlepoint, rand_pgnormalapprox

@enum PGSamplingMethod DEVROYE SADDLEPOINT NORMALAPPROX GAMMASUM DEVROYEPLUSGAMMASUM 

include("rand_pgdevroye.jl")
include("rand_pgsaddlepoint.jl")
include("rand_pgnormalapprox.jl")
include("rand_pggammasum.jl")
include("pghybrid.jl")

end
