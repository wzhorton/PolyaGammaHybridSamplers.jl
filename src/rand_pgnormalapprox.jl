"""
    rand_pgnormalapprox(b::Real, z::Real, [rng::AbstractRNG = Random.default_rng()])

Sample from a Polya-Gamma distribution using the normal approximation method.

# Arguments
- `b::Real`: the shape parameter
- `z::Real`: the exponential tilting parameter
- `rng::AbstractRNG`: random number generator object for `rand`

# Returns
- A sample from the Polya-Gamma distribution with shape parameter `b` and exponential tilting parameter `z`.

# Notes
- This method is an approximation that supports non-integer `b` and is very efficient for large `b`.
- Automatically selects this method when `b >= 170`.
"""
function rand_pgnormalapprox(b::Real, z::Real, rng::AbstractRNG)
    mu = mean(PolyaGammaHybridSampler(b, z))
    sigma = sqrt(var(PolyaGammaHybridSampler(b, z)))
    return rand(rng, Normal(mu, sigma))
end

function rand_pgnormalapprox(b::Real, z::Real)
    return rand_pgnormalapprox(b, z, Random.default_rng())
end