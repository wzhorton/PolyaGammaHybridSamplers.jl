"""
    PolyaGammaHybridSampler(b::Integer, z::Real, method::AbstractString)

A sampler for the Polya-Gamma distribution using a hybrid of the saddlepoint approximation, the Devroye method, 
and the normal approximation, each of which are discussed in Windle et al. (2014) -- see README.md for details.
The method used is determined by the the `method` parameter, which defaults to `"hybrid"` when not specified.

# Arguments
- `b::Integer`: The shape parameter of the Polya-Gamma distribution. Must be positive.
- `z::Real`: The exponential tilting parameter of the Polya-Gamma distribution.
- `method::PGSamplingMethod : An Enum object specifying the method used to sample from the Polya-Gamma distribution. 
    Must be one of `HYBRID`, `DEVROYE`, `SADDLEPOINT`, or `NORMALAPPROX`. Defaults to `HYBRID` when not specified.

# Returns
- A `PolyaGammaHybridSampler` object which can be sampled using `rand` or `rand!`.

# Examples
```julia
julia> using PolyaGammaHybridSamplers
julia> s = PolyaGammaHybridSampler(1, 1.0)
julia> rand(s)
```
"""
struct PolyaGammaHybridSampler{T <: Real, N <: Integer} <: Sampleable{Univariate, Continuous}
    b::N
    z::T
    method::PGSamplingMethod

    function PolyaGammaHybridSampler(b::N, z::T, method::PGSamplingMethod) where {T <: Real, N <: Integer}
        if b < zero(b) 
            throw(DomainError("b must be positive"))
        end
        new{T,N}(b, z, method)
    end
end


#----------------------------------#
# Define sampler objects
#----------------------------------#

# Define the outer constuctors 
function PolyaGammaHybridSampler(b::Integer, z::Real)
    PolyaGammaHybridSampler(b, z, HYBRID)
end

# Define mean and variance
pg_mean = function(b::Real, z::Real)
    if iszero(z)
        return b / 4.0
    else
        return b * inv(2.0*z) * tanh(z/2.0)
    end
end

pg_var = function(b::Real, z::Real)
    if iszero(z)
        return b * inv(24.0)
    else
        return b * inv(4.0*z^3) * (sinh(z)-z) * sech(z/2.0)^2
    end
end

Distributions.mean(s::PolyaGammaHybridSampler) = pg_mean(s.b, s.z)
Distributions.var(s::PolyaGammaHybridSampler) = pg_var(s.b, s.z)

#----------------------------------#
# Define rand hybrid methods
#----------------------------------#

function Base.rand(rng::AbstractRNG, s::PolyaGammaHybridSampler)
    if s.method == HYBRID
        return rand_pghybrid(s.b, s.z, rng)
    elseif s.method == DEVROYE
        return rand_pgdevroye(s.b, s.z, rng)
    elseif s.method == SADDLEPOINT
        return rand_pgsaddlepoint(s.b, s.z, rng)
    else # sampler.method = NORMALAPPROX
        return rand_pgnormalapprox(s.b, s.z, rng)
    end
end


"""
    rand_pghybrid(b::Integer, z::Real, rng::AbstractRNG)

Sample from a Polya-Gamma distribution using the hybrid selection method.

# Arguments
- `b::Integer`: the shape parameter
- `z::Real`: the exponential tilting parameter
- `rng::AbstractRNG`: random number generator object for `rand`

# Returns
- A sample from the Polya-Gamma distribution with shape parameter `b` and exponential tilting parameter `z`.

# Notes
- This method depends on the other rand_* methods in this package.
- Selection criteria: 
    - `b > 170` -> `rand_pgnormalapprox`
    - `b >= 13` -> `rand_pgsaddlepoint`
    - `b > 0` -> `rand_pgdevroye`
    - `b = 0` -> degenerate distribution at 0
- This function is not exported.
"""
function rand_pghybrid(b::Integer, z::Real, rng::AbstractRNG)
    if b > 170
        return rand_pgnormalapprox(b, z, rng)
    elseif b >= 13
        return rand_pgsaddlepoint(b, z, rng)
    elseif b > 0
        return rand_pgdevroye(b, z, rng)
    else # b = 0 -> degenerate distribution at 0
        return zero(b)
    end
end