"""
    rand_pgsaddlepoint(b::Real, z::Real, rng::AbstractRNG)

Sample from a Polya-Gamma distribution using the saddlepoint approximation method.

# Arguments
- `b::Real`: the shape parameter
- `z::Real`: the exponential tilting parameter
- `rng::AbstractRNG`: random number generator object for `rand`

# Returns
- A sample from the Polya-Gamma distribution with shape parameter `b` and exponential tilting parameter `z`.

# Notes
- This method is an approximation, but very efficient for large `b`.
- This method is recommended for `b > 13`, however no warning is given if `b` is too small.
- This method supports non-integer `b`.
"""

function rand_pgsaddlepoint(b::Real, z::Real, rng::AbstractRNG)
    if iszero(b) # b = 0 -> degenerate distribution at 0
        return zero(z)
    elseif b < 1
        throw(DomainError(b, "Saddlepoint approximation is only valid for b >= 1"))
    end

    return rand_Jstar_tilted(b, 0.5*abs(z), rng) * 0.25
end

"""
Note that we are sampling from the tilted J^*(b,z)

phi_z(x) = concave dual of K_z(t) = K(s) - s*x where s is plugged in to minimize the expression.
K(t) = cumulant generating function = log MGF = log E[exp(tX)]
M(t) = cosh^b(z/2) / cosh^b(sqrt((c^2/2 - t)/2) = cosh^b(z/2) * cosh^-b(1/2 * sqrt(c^2 - 2t))
K'(t) = (b/2)*tanh(1/2 * sqrt(c^2 - 2t)) / sqrt(c^2 - 2t)

m = mode of phi_z(x) = argmax_x phi_z(x) = tanh(z)/z

sp_n(x|z) = saddle point approximation of f.
    = sqrt(n/2pi) K''(t(x))^(-1/2) * exp(n*phi_z(x))
where t(x) = argmax_t K(t) - t*x
    t(x) = 1/2 * sqrt(c^2 - 2x)
    K''(t) = (b/4) * (c^2 - 2t) * sech^2(1/2 * sqrt(c^2 - 2t)) / (c^2 - 2t)^(3/2)
    K''(t(x)) = (b/4) * sech^2(1/2 * sqrt(c^2 - 2x)) / (c^2 - 2x)
    sp_n(x|z) = sqrt(n/2pi) * (b/4)^(1/2) * sech(1/2 * sqrt(c^2 - 2x)) / (c^2 - 2x)^(3/4) * exp(n*phi_z(x))

eta(x) = phi(x) - delta(x)
delta(x) = 1/2*(1/xc + 1/x) if x <= xc, otherwise log(x) - log(xc)

k(x|n,z) = sqrt(n/2pi) * cases:
    x < xc: 1/sqrt(alpha_l) * exp(n/(2xc)) * x^(-3/2) * exp(-n/(2x) + nLl(x|z))
    x > xc: 1/sqrt(alpha_r) * xc^n * x^(n-1) * exp(nLr(x|z))

the left side is inverse Gaussian while the right is gamma.
=> k_left is inverse Gaussian with mu = 1/sqrt(pl) and lambda = n,
    pl = -2Ll'(x)
   k_right is gamma with shape = n, rate = npr
    pr = -Lr'(x)

MAIN ALGO:
xl = m, xc = 1.1 * m, xr = 1.2 * m
calculate the tangent lines of eta at xl and xr: call them Ll(x|z) and Lr(x|z).
proposal g is propto k which is mixture of 2 parts: inverse Gaussian and gamma.

Draw X from g,
draw U from Unif(0, k(X))
if U < sp_n(X|z), accept X, otherwise reject X.

    ---

a lot of what I wrote above is slightly inaccurate. The Jstar tilted distribution is what we are sampling from
and it has a MGF of cosh^b(z)*cosh^-b(sqrt(z^2 - 2t)).
"""

# The following code is almost a direct translation of the C++ code from the original BayesLogit package,
# which is available at https://github.com/jwindle/BayesLogit/blob/master/src/
# To help with organization, helper functions are organized by the file they originally appeared in.

#---------------------------------#
# Main function for Jstar tilted
# PolyaGammaApproxSP.cpp
#---------------------------------#

function rand_Jstar_tilted(b::Real, z::Real, rng::AbstractRNG)
    max_iter = 1000
    xl = y_eval(-z^2)
    md = xl * 1.1
    xr = xl * 1.2

    vmd = v_eval(md)
    K2md = 0.0
    if abs(vmd) >= 1e-6
        K2md = md^2 + (1-md) *inv(vmd)
    else
        K2md = md^2 - 1/3 - (2/15) * vmd
    end

    m2 = md^2
    al = m2 * md * inv(K2md)
    ar = m2 * inv(K2md)

    rl, il = tangent_to_eta(xl, z, md)
    rr, ir = tangent_to_eta(xr, z, md)

    rl = -rl
    rr = -rr

    lcn = 0.5 * log(0.5 * b *inv(π))
    rt2rl = sqrt(2 * rl)

    wl = exp(0.5 * log(al) - b * rt2rl + b * il + 0.5 * b * inv(md)) *
         p_igauss(md, inv(rt2rl), b)

    wr = exp(0.5 * log(ar) + lcn + (-b * log(b * rr) + b * ir - b * log(md) + lgamma(b))) *
         (1.0 - cdf(Gamma(b, inv(b*rr)), md))

    wt = wl + wr
    pl = wl * inv(wt)

    go = true
    iter = 0
    X = 2.0
    F = 0.0

    while go && iter < max_iter
        iter += 1

        if rand(rng) < pl
            X = rtigauss(inv(rt2rl), b, md, rng)
            phi_ev = b * (il - rl * X) + 0.5 * b * ((1-inv(X)) - (1-inv(md)))
            F = exp(0.5 * log(al) + lcn - 1.5 * log(X) + phi_ev)
        else
            X = ltgamma(b, b*rr, md, rng)
            phi_ev = b * (ir - rr * X) + b * (log(X) - log(md))
            F = exp(0.5 * log(ar) + lcn + phi_ev) * inv(X)
        end

        spa = sp_approx(X, b, z)

        if F * rand(rng) < spa
            go = false
        end
    end

    return b * X
end

#--------------------------------#
# Helper Functions: InvertY.cpp
#--------------------------------#

const YGRID = [0.0625,0.06698584,0.07179365,0.07694653,0.08246924,
0.08838835,0.09473229,0.1015315,0.1088188,0.1166291,
0.125,0.1339717,0.1435873,0.1538931,0.1649385,
0.1767767,0.1894646,0.2030631,0.2176376,0.2332582,
0.25,0.2679434,0.2871746,0.3077861,0.329877,
0.3535534,0.3789291,0.4061262,0.4352753,0.4665165,
0.5,0.5358867,0.5743492,0.6155722,0.659754,
0.7071068,0.7578583,0.8122524,0.8705506,0.933033,
1,1.071773,1.148698,1.231144,1.319508,
1.414214,1.515717,1.624505,1.741101,1.866066,
2,2.143547,2.297397,2.462289,2.639016,
2.828427,3.031433,3.24901,3.482202,3.732132,
4,4.287094,4.594793,4.924578,5.278032,
5.656854,6.062866,6.498019,6.964405,7.464264,
8,8.574188,9.189587,9.849155,10.55606,
11.31371,12.12573,12.99604,13.92881,14.92853,
16]

const VGRID = [-256,-222.8609,-194.0117,-168.897,-147.0334,
-128,-111.4305,-97.00586,-84.4485,-73.51668,
-63.99997,-55.71516,-48.50276,-42.22387,-36.75755,
-31.99844,-27.85472,-24.24634,-21.10349,-18.36524,
-15.97843,-13.89663,-12.07937,-10.49137,-9.101928,
-7.884369,-6.815582,-5.875571,-5.047078,-4.315237,
-3.667256,-3.092143,-2.580459,-2.124095,-1.716085,
-1.350442,-1.022007,-0.7263359,-0.4595871,-0.2184366,
0,0.1982309,0.3784427,0.5425468,0.6922181,
0.828928,0.953973,1.068498,1.173516,1.269928,
1.358533,1.440046,1.515105,1.584282,1.64809,
1.706991,1.761401,1.811697,1.858218,1.901274,
1.941143,1.978081,2.012318,2.044068,2.073521,
2.100856,2.126234,2.149802,2.171696,2.192042,
2.210954,2.228537,2.244889,2.260099,2.274249,
2.287418,2.299673,2.311082,2.321703,2.331593,
2.340804]

function y_eval(v)
    tol = 1e-6
    y = 0.0
    r = sqrt(abs(v))
    if v > tol
        y = tan(r) / r
    elseif v < -tol
        y = tanh(r) / r
    else
        y = 1 + (1/3) * v + (2/15) * v^2 + (17/315) * v^3
    end
    return y
end

function ydy_eval(v)
    tol = 1e-6
    y = y_eval(v)
    if abs(v) >= tol 
        yp = y
        dyp = 0.5 * (y^2 + (1-y) * inv(v))
    else
        yp = y
        dyp = 0.5 * (y^2 - (1/3) - (2/15) * v)
    end
    return yp, dyp
end

function f_eval(v, y)
    return y_eval(v) - y
end

function fdf_eval(v, y)
    fp, dfp = ydy_eval(v)
    fp -= y
    return fp, dfp
end

function df_eval(v)
    f, df = ydy_eval(v)
    return df
end

function v_eval(y)
    ylower = YGRID[1]
    yupper = YGRID[end]
    tol = 1e-9
    max_iter = 1000

    if y < ylower
        return -inv(y^2)
    elseif y > yupper
        return atan(0.5 * π * y)^2
    elseif isone(y)
        return 0
    end

    id = 10*(log(y) / log(2.0) + 4.0)
    idlow = Int(floor(id))
    idhigh = Int(ceil(id))
    vl = VGRID[idlow]
    vh = VGRID[idhigh]

    iter = 0
    diff = tol + 1.0
    vnew = vl
    vold = vl

    while diff > tol && iter < max_iter
        iter += 1
        f0, f1 = fdf_eval(vold, y)
        vnew = vold - f0 / f1
        vnew = min(vnew, vh)
        vnew = max(vnew, vl)
        diff = abs(vnew - vold)
        vold = vnew
    end

    if iter >= max_iter
        @error "v_eval: reached max_iter: $iter"
    end

    return vnew
end


#---------------------------------#
# Helper functions: PolyaGammaApproxSP.cpp
#---------------------------------#

function cos_rt(v)
    r = sqrt(abs(v))
    if v >= 0
        y = cos(r)
    else
        y = cosh(r)
    end
    return y
end

function delta_func(x, mid)
    if x >= mid
        delta_val = log(x) - log(mid)
        delta_der = inv(x)
    else
        delta_val = 0.5 * (1 - inv(x)) - 0.5 * (1 - inv(mid))
        delta_der = 0.5 * inv(x^2)
    end
    return delta_val, delta_der
end

function phi_func(x, z)
    v = v_eval(x)
    u = 0.5 * v
    t = u + 0.5 * z^2

    phi_val = log(cosh(abs(z))) - log(cos_rt(v)) - t * x
    phi_der = -t

    return phi_val, phi_der
end

function tangent_to_eta(x, z, mid)
    phi_val, phi_der = phi_func(x, z)
    delta_val, delta_der = delta_func(x, mid)

    eta_val = phi_val - delta_val
    eta_der = phi_der - delta_der

    slope = eta_der
    intercept = eta_val - eta_der * x

    return slope, intercept
end

function sp_approx(x, n, z)
    v = v_eval(x)
    u = 0.5 * v
    t = u + 0.5 * z^2

    phi = log(cosh(z)) - log(cos_rt(v)) - t * x

    K2 = 0.0
    if abs(v) >= 1e-6
        K2 = x^2 + (1-x) * inv(v)
    else
        K2 = x^2 - 1/3 - v * (2/15)
    end
    log_spa = 0.5 * log(0.5 * n / pi) - 0.5 * log(K2) + n * phi
    return exp(log_spa)
end

function rtigauss(mu, lambda, trunc, rng::AbstractRNG)
    X = trunc + 1.0
    if trunc < mu
        alpha = 0.0
        while rand() > alpha
            X = rtinvchi2(lambda, trunc, rng)
            alpha = exp(-0.5 * lambda * inv(mu^2) * X)
        end
    else
        while X > trunc
            X = igauss(mu, lambda, rng)
        end
    end
    return X
end


#---------------------------------#
# Helper functions: inverse_gaussian.cpp
#---------------------------------#

function igauss(mu, lambda, rng::AbstractRNG)
    mu2 = mu^2
    Y = rand(rng, Normal())^2
    W = mu + 0.5 * mu2 * Y * inv(lambda)
    X = W - sqrt(W^2 - mu2)
    if rand(rng) > mu / (mu + X)
        X = mu2 * inv(X)
    end
    return X
end

function p_igauss(x, mu, lambda)
    z = inv(mu)
    b = sqrt(lambda * inv(x)) * (x*z - 1)
    a = sqrt(lambda * inv(x)) * (x*z + 1) * -1.0
    y = cdf(Normal(),b) + exp(2 * lambda * z + logcdf(Normal(),a))
    return y
end


#---------------------------------#
# Helper functions: truncated_norm.cpp
#---------------------------------#

function rtinvchi2(scale, trunc, rng::AbstractRNG)
    R = trunc / scale
    E = rand(rng, truncated(Normal(); lower = inv(sqrt(R))))
    X = scale * inv(E^2)
    return X
end


#---------------------------------#
# Helper functions: truncated_gamma.cpp
#---------------------------------#

function ltgamma(shape, rate, trunc, rng::AbstractRNG)
    a = shape
    b = rate * trunc

    if shape < one(shape)
        @error "ltgamma: shape < 1"
    end

    if isone(shape)
        return rand(rng, Exponential()) / rate + trunc
    end

    d1 = b - a
    d3 = a - 1
    c0 = 0.5 * (d1 + sqrt(d1^2 + 4b)) / b

    x = 0.0
    accept = false

    while !accept
        x = b + rand(rng, Exponential()) / c0

        l_rho = d3 * log(x) - x * (1 - c0)
        l_M   = d3 * log(d3 / (1 - c0)) - d3

        accept = log(rand(rng)) <= (l_rho - l_M)
    end

    return trunc * (x / b)
end