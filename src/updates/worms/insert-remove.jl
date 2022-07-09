"""
    Insert(s::System, slices::Int64, <keyword arguments>)

Initiates the Insert struct where the variables for the Insert update is stored. The Insert update is one of the updates of the Worms Algorithm [1]. It inserts worms when the simulation is in the Z-sector and thereby transform the configuration the the G-sector. It is completempary to the Remove update.

# Keyword arguments
- `minslices = 1`: The minimum slices the Insert-function can add into the simulation. Must be larger or equel one.
- `maxslices = s.M*s.N`: The maximum slices the Insert-function can can add into the simulation. These variable can become arbitrary long. The default is set to the number of time sliced of a worldline times the number of particle.
- `minacc=0.4`: The lower bound accaptance used to change the number of slices inserted by the Insert-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices inserted by the Insert-function.

[1] 10.1103/PhysRevE.74.036701
"""

struct Insert <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function Insert(d::Threshold; adj=10, range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::Insert)(s::System)::Bool
    acc = false
    # only if config in z-sector and worms algorithm is turned on
    if s.worms > 0 || s.worm_algorithm == false
        return acc
    end

    # choose the add wich time slice `j₀` and how much beads `m` are added above
    j₀ = rand(1:s.M)
    m = rand(1:u.var.m) # (m+1)-bead long worm; m links
    jₘ = j₀ + m # last bead

    # generate the new beads with the levy construction
    # the new beads must be further than the 2d scattering length from any other bead (hard sphere)
    r′ = zeros(Float64, m + 1, s.dim) # list of new positions
    # generate first bead randomly and ensure it is not in a hardsphere of another bead
    r′[begin, :] = uniform_sample_hardsphere(s, j₀)
    if all(r′[begin, :] .== 0.0)
        queue!(u.counter_var, acc)
        return acc
    end

    # generate last bead close enough to the first, i.e., √(m * s.τ), so we have a high accaptance rate
    # again the new bead cannot lie another hardsphere
    r′[end, :] = gaussian_sample_hardsphere(s, r′[begin, :], jₘ, m)
    if all(r′[end, :] .== 0.0)
        queue!(u.counter_var, acc)
        return acc
    end

    # the levy construction computes remaining beads
    bool = hardspherelevy!(r′, s, j₀)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # copmute the relevant links and the action
    w_updated = 0.0 # the action of the updated configuration
    V′ = zeros(Float64, m) # new links
    for j in 1:m
        # for every timeslice add potential contribution
        V′[j] = s.lnV(r′[j, :], r′[j+1, :])

        # For all NeirestNeighbours compute the interaction contribution
        interaction_action!(w_updated, s, r′, j, j₀)
    end
    # updated weight; there is no initial weight as only beads are added, not changed
    w_updated += sum(V′) # updated weight

    preamble = get_C(s.C) * u.var.m * s.vol * s.M
    # metropolis question
    if metropolis(preamble * exp(w_updated + s.μ * m * s.τ))
        acc = true

        nw = ((jₘ - 1) ÷ s.M) + 1
        jj, ll = 1, 1
        for i in 1:nw
            r = Matrix{Union{Missing,Float64}}(missing, s.M, s.dim)
            V = Vector{Union{Missing,Float64}}(missing, s.M)
            bins = Vector{Union{Missing,Int64}}(missing, s.M)

            # add computed positions and links
            for j in jj:(m+1)
                rⱼ = @view r′[j, :]
                r[mod1(j₀ + j - 1, s.M), :] = rⱼ
                bins[mod1(j₀ + j - 1, s.M)] = bin(rⱼ, s.nbins, s.L)

                if mod1(j₀ + j - 1, s.M) == s.M
                    jj = j + 1
                    break
                end
            end
            for j in ll:m
                V[mod1(j₀ + j - 1, s.M)] = V′[j]
                if mod1(j₀ + j - 1, s.M) == s.M
                    ll = j + 1
                    break
                end
            end

            # sidecases get indexing right
            if i == 1
                worm = true
                if nw != 1
                    s.worms = 2
                    tail = j₀
                    head = 0
                    next = nw == 2 ? s.N + 2 : s.N + 1
                    add_nn!(s.nn,  s.N + nw - 1, bins)
                else
                    s.worms = 1
                    tail = j₀
                    head = jₘ
                    next = 0
                    add_nn!(s.nn, s.N + 1, bins)
                end
            elseif i == nw
                worm = true
                tail = 0
                head = mod1(jₘ, s.M)
                next = 0
                add_nn!(s.nn,  s.N + nw, bins)
            else
                worm = false
                next = i + 1 == nw ? length(s.world) + 2 : s.N + i
                add_nn!(s.nn, s.N + i -1, bins)
            end

            # insert worms or particle
            if worm
                push!(s.world, Worm(r, V, bins, tail, head, next))
            else
                insert!(s.world, length(s.world), Particle(r, V, bins, next))
            end
        end
        s.Nz = s.N
        s.N = max(s.N + nw - 2, s.N)
    end
    queue!(u.counter_var, acc)
    return acc
end

"""
    Remove(s::System, slices::Int64, <keyword arguments>)

Initiates the Remove struct where the variables for the Remove update is stored. The Remove update is one of the updates of the Worms Algorithm [1]. It removes worms when the simulation is in the G-sector and thereby make a z-sector configuration. It is completempary to the Insert update. The second argumant `slices` is used as a threshold. If the worm chosen has more bead as this threshold the update is rejected.

# Keyword arguments
- `minacc=0.4`: The lower bound accaptance used to change the number of slices inserted by the Remove-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices inserted by the Remove-function.

[1] 10.1103/PhysRevE.74.036701
"""

struct Remove <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function Remove(
        d::Threshold;
        adj=10,
        range=100_000
    )
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::Remove)(s::System)::Bool
    acc = false
    # only if config in G-sector and the Worms algorithm is turned on
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end

    # find the tail of the Worm
    ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))

    # compute the length of the polymer of the worm
    Npol, pol = subcycle(s.world, ntail)
    m = lengthpolymer(s.world, pol, :link) # number of links

    # If the number of particle are changed when we work in the canonical ensemble, reject
    if !s.gc && s.N - max(Npol - 2, 0) != s.Nz
        return acc
    end

    # if the length of the polymer is larger than the number of slices the move can remove, we reject the move
    # the number of slices the move can remove is changed according the accaptance ratio
    if m > u.var.m
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the inital action; there is no updated action because beads are only removed
    w_initial = 0.0 # initial weight
    for n in pol
        # add all potential contributions to the action
        w_initial += sum(skipmissing(s.world[n].V))

        for j in eachindex(skipmissing(view(s.world[n].r, :, 1)))
            # add all nn interaction to the action
            interaction_action!(w_initial, s, _ -> n, j)
        end
    end

    preamble = (1 / (get_C(s.C)* u.var.m * s.vol * s.M))
    # metropolis question
    if metropolis(preamble * exp(-w_initial - s.μ * m * s.τ))
        acc = true

        rm_nn!(s.nn, pol, s.world)
        # if in polymer, than fix the indexing in s.world and delete
        fixindexing!(s.world, s.nn, pol)
        deleteat!(s.world, sort(pol))
        s.N -= max(Npol - 2, 0)
        s.worms = 0
    end
    queue!(u.counter_var, acc)
    return acc
end


