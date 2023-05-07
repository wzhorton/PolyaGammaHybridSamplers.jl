using Test
using Random
using Distributions
using SpecialFunctions
using PolyaGammaHybridSamplers
import PolyaGammaHybridSamplers: PGSamplingMethod

include("pghybrid_tests.jl")
include("rand_pgnormalapprox_tests.jl")
include("rand_pgdevroye_tests.jl")
include("rand_pgsaddlepoint_tests.jl")