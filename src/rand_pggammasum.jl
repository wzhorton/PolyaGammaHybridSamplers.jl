"""
    rand_pggammasum(b::Real, z::Real, rng::AbstractRNG)

Sample from a Polya-Gamma distribution using the truncated sum of gammas representation.

# Arguments
- `b::Real`: the shape parameter
- `z::Real`: the exponential tilting parameter
- `rng::AbstractRNG`: random number generator object for `rand`

# Returns
- A sample from the Polya-Gamma distribution with shape parameter `b` and exponential tilting parameter `z`.

# Notes
- This method is an approximation, meant only for b < 1
- No warning is given if `b` is too large.
- This method supports non-integer `b`.
- This function is not exported.
"""
function rand_pgpggammasum(b::Real, z::Real, rng::AbstractRNG)
    if iszero(b) # b = 0 -> degenerate distribution at 0
        return zero(z)
    end
    trunc_level = 200 # see paper
    total = zero(z)
    for k in 1:trunc_level
        total += rand(rng,Gamma(b, 1.0)) / (z^2 + (k - 0.5)^2 * 4 * π^2)
    end
    return 2 * total
end