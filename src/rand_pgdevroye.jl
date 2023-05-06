"""
    rand_pgdevroye(b::Integer, z::Real, rng::AbstractRNG)

Sample from a Polya-Gamma distribution using the Devroye method.

# Arguments
- `b::Integer`: the shape parameter
- `z::Real`: the exponential tilting parameter
- `rng::AbstractRNG`: random number generator object for `rand`

# Returns
- A sample from the Polya-Gamma distribution with shape parameter `b` and exponential tilting parameter `z`.

# Notes
- This method is exact, but is increasingly slower as `b` increases.
"""

function rand_pgdevroye(b::Integer, z::Real, rng::AbstractRNG)
    if iszero(b) # b = 0 -> degenerate distribution at 0
        return zero(z)
    end 

    draw  = zero(z)
    for _ in Base.OneTo(b)
        draw += rand_Jstar(0.5*abs(z), rng) * 0.25
    end
    return draw
end


#----------------------------------#
# Define auxiliary functions 
#----------------------------------#

# Define the J* sampler function (See Algorithm 1, Supplementary Materials of Polson et al. 2013)
function rand_Jstar(z::Real, rng::AbstractRNG)
    t = 0.64 # paper recommends this constant
    K = 0.125*π^2 + 0.5*z^2
    p = 0.5*π*inv(K) * exp(-t*K)
    q = 2*exp(-z) * cdf(InverseGaussian(inv(z), 1.0), t)
    while true
        # Generate X
        if rand(rng, Uniform(0.0, 1.0)) < p/(p+q) #U ~ Uniform(0,1)
            # Draw exponential(K) truncated above t
            X = t + inv(K)*rand(rng, Exponential())
        else
            # Draw inverseGaussian(1/z, 1) truncated below t
            # WARNING: This implementation is unstable. For z ≈ 0, an infinite loop occurs
            # due to the InverseGaussian parameter being infinite. For large z (>400 ish),
            # the loop also never terminates, probably due to underflow in a_coefs().
            # z + 0.0001 is a hack to address the case where z = 0. In the future,
            # Algorithms 2 and 3 of Polson et al. 2013 should be implemented to fix this.
            X = rand(rng, truncated(InverseGaussian(inv(z+0.0001), 1.0); upper = t))
        end
        # Accumulate a(X) to S
        S = a_coefs(0, X, t)
        Y = S * rand(rng, Uniform(0.0, 1.0)) #V ~ Uniform(0,1)
        n = Int(0)
        while true
            n += 1
            if isodd(n)
                S -= a_coefs(n, X, t)
                if Y < S
                    return X # accept X and exit
                end
            else
                S += a_coefs(n, X, t)
                if Y > S
                    break # reject X and start over
                end
            end
        end
    end
end

# Define the a_coefs function (See Equations (15, 16), Polson et al. 2013)
function a_coefs(n::Integer, x::Real, t::Real)
    common_terms = π * (n+0.5) * exp(-(n+0.5)^2)
    if x <= t
        return common_terms * exp(-2 * inv(x))
    else 
        return common_terms * exp(-x * 0.5 * π^2)
    end
end