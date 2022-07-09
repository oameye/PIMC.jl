# As in every update the propagators from  R_1 to R_2 is greatly reduces in the hastings question we do not write a function which represent the small time propagator completly. Instead we only write function for its componenents

Coord = Union{Vector{Float64}, SubArray{Float64, 1}}

# Propagator helper function for calculating distance over periodic boundary conditions
function distance(x1::Float64, x2::Float64, L::Float64)::Float64
    dx = abs(x1 - x2)
    min((2 * L) - dx, dx) # benched
end

# function distance(x1::Float64, x2::Float64, L::Float64)::Float64
#     return PeriodicEuclidean(2*L)(x1 + L, x2 + L)
# end

# Only completely used in constructor of the struct system, afterward a curry function
function lnK(r1::Coord, r2::Coord, τ::Float64, λ::Float64, L::Float64)::Float64
    dr = distance.(r1, r2, L) # vector
    - dot(dr, dr) / (4*λ*τ)
end
function prop_0(r1::Coord, r2::Coord, τ::Float64, λ::Float64, L::Float64, dim::Int64)::Float64
    dr = distance.(r1, r2, L) # vector
    # return  (1/√(4*π*λ*τ))^dim * exp(- dot(dr, dr) / (4*λ*τ))
    return exp(- dot(dr, dr) / (4*λ*τ))
end

function lnV(r1::Coord, r2::Coord, τ::Float64, V::Function)::Float64
    -0.5 * τ * (V(r1) + V(r2))
end

function teleport(x::Float64, L::Float64)::Float64
    (x + L) - floor(x / (2 * L) + 0.5) * (2 * L) - L
end

# the integral terms of prop_rel need only two scalars; to not generate a 4 dimensional interpolation function we only interpolate the integral terms
function prop_rel_interpolate_terms(L::Float64, g0::Float64, τ::Float64)
    Δ = 600

    # Array of norm of positionvectors r
    A_r1 = range(1e-20, L, Δ) # boundary of periodic boundary condition (PBC)
    A_r2 = range(1e-20, L, Δ) # r will always be smaller than L because of PBC
    # => the integrals blow up for r=0 and cannot be evaluated

    # helpfunction
    tk(k) = 1 / ((2 / pi) * (Base.MathConstants.eulergamma + log(k / 2)) - (4 / g0))

    # Integral term (Gautier et al. supplementary information S28)
    term1(r1_norm, r2_norm) = (1 / (2 * pi)) * quadgk(k ->
    k * exp(-τ * k^2  ) * (tk(k)^2 / (1 + tk(k)^2)) *
    besselj0(k * r1_norm) * besselj0(k * r2_norm)
    , 0, Inf, rtol = 1e-11, order = 15, maxevals=10^7)[1]

    term2(r1_norm, r2_norm) = (1 / (2 * pi)) * quadgk(k ->
    k * exp(-τ * k^2  ) * (tk(k) / (1 + tk(k)^2)) *
    (besselj0(k * r1_norm) * bessely0(k * r2_norm) + besselj0(k * r2_norm) * bessely0(k * r1_norm))
    , 0, Inf, rtol = 1e-11, order = 15, maxevals=10^7)[1]

    term3(r1_norm, r2_norm) = (1 / (2 * pi)) * quadgk(k ->
    k * exp(-τ * k^2) * (tk(k)^2 / (1 + tk(k)^2)) *
    bessely0(k * r1_norm) * bessely0(k * r2_norm)
    , 0, Inf, rtol = 1e-11, order = 15, maxevals=10^7)[1]

    int_terms(r1_norm, r2_norm) = - term1(r1_norm, r2_norm) - term2(r1_norm, r2_norm) + term3(r1_norm, r2_norm)

    A = [int_terms(r1_norm, r2_norm) for r1_norm in A_r1, r2_norm in A_r2]

    #scaled BSpline interpolation
    sitp = scale(interpolate(A, BSpline(Linear())), A_r1, A_r2)

    return sitp
end

# the relative free space propagator; note 'relative' hence m → μ=m/2
function prop_rel0(r1_rel::Vector, r2_rel::Vector, τ::Float64)
    # Computes the free relative small time propagator of Gautier et al.
    return exp(-dot(r1_rel - r2_rel, r1_rel - r2_rel) / (4 * τ))/(4 * π * τ)
end

# the interaction propagator
function build_prop_int(L::Float64, g0::Float64, τ::Float64)
    terms = prop_rel_interpolate_terms(L, g0, τ)

    function prop_int(r1_rel::Vector, r2_rel::Vector, τ::Float64)
        r1_norm = norm(r1_rel)
        r2_norm = norm(r2_rel)
        return 1+terms(r1_norm, r2_norm)/prop_rel0(r1_rel, r2_rel, τ)
    end

    return prop_int
end

