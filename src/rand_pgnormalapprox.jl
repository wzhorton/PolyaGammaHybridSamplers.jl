"""
    rand_pgnormalapprox(b::Real, z::Real, rng::AbstractRNG)

Sample from a Polya-Gamma distribution using the normal approximation method.

# Arguments
- `b::Real`: the shape parameter
- `z::Real`: the exponential tilting parameter
- `rng::AbstractRNG`: random number generator object for `rand`

# Returns
- A sample from the Polya-Gamma distribution with shape parameter `b` and exponential tilting parameter `z`.

# Notes
- This method is an approximation, but very efficient for large `b`.
- This method is recommended for `b > 170`, however no warning is given if `b` is too small.
- This method supports non-integer `b`.
- This function is not exported.
"""
function rand_pgnormalapprox(b::Real, z::Real, rng::AbstractRNG)
    mu = pg_mean(b, z)
    sigma = sqrt(pg_var(b, z))
    return rand(rng, Normal(mu, sigma))
end