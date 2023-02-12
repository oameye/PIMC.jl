mutable struct Particle <: Worldline
    r::Matrix{Float64}
    V::Vector{Float64}
    bins::Vector{Int64}
    next::Int64
end
# TODO Matrix of dictionairies
NearestNeighbors = Vector{Vector{Vector{Int64}}}

function determine_nnrange(propint::Function, τ::Float64, a::Float64, b::Float64)::Float64
    onevarpropint(r::Float64) = propint([r], [r], τ) - 0.999
    opt = optimize(x -> onevarpropint(first(x)), [a])
    r_min = Optim.minimizer(opt)
    return find_zero(onevarpropint, (r_min[1], b), Bisection())
end

function init_int(interactions::Bool, propint::Function, τ::Float64, L::Float64, rₐ::Float64, λ::Float64, μ::Float64)::Tuple{Float64,Function}
    if !interactions
        lnU = function (r1_rel::Vector{Float64}, r2_rel::Vector{Float64}, τ::Float64)
            0.0
        end
        if rₐ == 0.0
            rₐ = L / 4 # length of box per dimension is 2*L
        end
    else
        lnU = function (r1_rel::Vector{Float64}, r2_rel::Vector{Float64}, τ::Float64)
            propint(r1_rel, r2_rel, τ) < 0.0 ? -μ : log(propint(r1_rel, r2_rel, τ))
        end
        if rₐ == 0.0
            rₐ = determine_nnrange(propint, τ, 1e-20, L)
        end
    end
    return (rₐ, lnU)
end

function init_world(potential::Function, N::Int64, M::Int64, dim::Int64, L::Float64, τ::Float64, a::Float64, λ::Float64)::Vector{Worldline}
    ps = Vector{Worldline}(undef, N)
    for n in eachindex(ps)
        V = Vector{Float64}(undef, M)
        r = Matrix{Float64}(undef, M, dim)

        ctr = 0
        pass = n == 1 ? false : true
        while pass
            pass = false
            ctr += 1
            if ctr > 10_000
                @error "creating world of $(N) particles with hardspheres of length $(a) in the volume $((2*L)^dim) failed!"
                return Vector{Worldline}()
            end

            r[begin, :] = 2 * L .* (rand(Float64, dim) .- 0.5)
            r[end, :] = r[begin, :]
            levy!(r, τ, L, λ)
            for m in eachindex(V)
                for i in 1:(n-1)
                    d = norm(distance.(ps[i].r[m, :], r[m, :], L))
                    if d < a
                        levy!(r, τ, L, λ)
                        pass = true
                    end
                end
            end
        end
        if n == 1
            r[begin, :] = 2 * L .* (rand(Float64, dim) .- 0.5)
            r[end, :] = r[begin, :]
            levy!(r, τ, L, λ)
        end

        # Pre-calc propagator for each time slice
        for m in eachindex(V)
            V[m] = lnV(r[m, :], r[mod1(m + 1, M), :], τ, potential)
        end
        ps[n] = Particle(r, V, zeros(M), n)
    end
    return ps
end

function init_nn(world::Vector{Worldline}, M::Int64, dim::Int64, L::Float64, rₐ::Float64)::Tuple{NearestNeighbors, Vector{Vector{Int64}}, Int64}
    nbins = round(Int64, (2 * L) / rₐ, RoundDown)
    # bin = 2*L/nbins # bin size
    nn = NearestNeighbors(undef, M)
    for m in eachindex(nn)
        bins = [Int64[] for _ in 1:nbins^dim]
        psort!(world, bins, nbins, m, L)
        nn[m] = bins
    end
    nbs = Vector{Int64}[bin_neighbors(i, nbins, dim) for i in 1:nbins^dim]
    return nn, nbs, nbins
end

mutable struct System
    dim::Int64 # dimensions
    M::Int64 # Number of time slices
    N::Int64 # Number of particles
    Ninit::Int64
    μ::Float64 # Chemical potential; units of Eᵣ/kᵦ
    λ::Float64 # units

    L::Float64 # radius of symmteric box around origen
    vol::Float64

    β::Float64 # Inverse temperature
    τ::Float64 # Time slice

    # the typical structure of world is: PPPTH (N=3, worms=2)
    world::Vector{Worldline}

    nn::NearestNeighbors # nearest noughbour grid
    nbs::Vector{Vector{Int64}} # neigbouring bins of bins
    nbins::Int64 # number of bins in nearest noughbour grid

    V::Function # potential function
    dV::Function # derivative of the potential
    lnV::Function # Curry the potential propagator function
    lnK::Function # Curry the kinetic propagator function

    lnU::Function # Curry the interaction propagator function
    a::Float64 # The cut-off hardsphere parameter (2D scattering length)
    ctr::Int64 # number of tries to sample beads outside hardsphere

    measure_scheme::Symbol # symbol wich indicates which measure scheme to use
    N_MC::Dict{Int64, Int64} # dictionairies of the number of measurement taken
    Nctr::Dict{Int64, Int64} # dictionairy of counter which determines if we take a measurement
    Ncycle::Int64 # after how much canonical configuration we take a measurement

    # Initialization system struct with a constructor
    function System(
        potential::Function;
        dV::Function = zero,
        dim::Int64 = 2,
        M::Int64 = 100,
        N::Int64 = 2,
        μ::Float64 = 0.0,
        L::Float64 = 4.0,
        T::Float64 = 1.0,
        λ::Float64 = 1.0,
        interactions::Bool = false,
        propint::Function = _ -> 0.0,
        g::Float64 = 0.0,
        rₐ::Float64 = 0.0,
        length_measurement_cycle::Int64 = 10,
        measure_scheme::Symbol = :c
    )

        β = 1.0 / T
        τ = β / M
        vol = (2 * L)^dim

        a = interactions ? exp(-2 * pi / g) : 0.0

        # Random initial positions for each particle
        world = init_world(potential, N, M, dim, L, τ, a, λ)

        rₐ, lnU = init_int(interactions, propint, τ, L, rₐ, λ, μ)

        nn, nbs, nbins = init_nn(world, M, dim, L, rₐ)

        new(dim, M, N, N, μ, λ, L, vol, β, τ, world, nn, nbs, nbins,
            potential, dV,
            (x1, x2) -> lnV(x1, x2, τ, potential),
            (x1, x2, τ) -> lnK(x1, x2, λ, τ, L),
            (x1, x2) -> lnU(x1, x2, τ),
            a, 10_000,
            measure_scheme, Dict{Int64, Int64}([(N, 0)]), Dict{Int64, Int64}([(N, 0)]), length_measurement_cycle)
    end
end
