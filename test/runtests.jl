using Test
using Random
using Distributions
using SpecialFunctions
using PolyaGammaHybridSamplers
import PolyaGammaHybridSamplers: PGSamplingMethod

include("PolyaGammaHybridSamplers_test.jl")
include("rand_pgnormalapprox_tests.jl")
include("rand_pgdevroye_tests.jl")