"""
    OpenGC(s::System, slices::Int64, <keyword arguments>)

Initiates the Open struct where the variables for the Open update is stored. The Open update is one of the updates of the Worms Algorithm [1]. It opens particles with a certain amount of beads(`slices`) when the simulation is in the Z-sector and thereby transform the configuration the the G-sector. It is completempary to the Close update.

# Keyword arguments
- `maxslices = s.M-1`: The maximum slices the Open-function can Open into the simulation.
- `minacc=0.4`: The lower bound accaptance used to change the number of slices opened by the open-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices opened by the Open-function.

[1] 10.1103/PhysRevE.74.036701
"""

struct OpenGC <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function OpenGC(d::Threshold; adj=10, range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::OpenGC)(s::System)::Bool
    acc = false
    # only if config in Z-sector
    if s.worms > 0 || s.worm_algorithm == false || s.N==0
        return acc
    end

    if s.N == 0
        error("N is equal to zeros but s.worms == $(s.worms). Weird?!")
    end

    # choose a bead random
    j₀, n = rand(1:s.M), rand(1:s.N)

    # determine the polymer
    Npol, pol = subcycle(s.world, n)
    fpcycle(j::Int64)::Int64 = pcycle(j, pol, Npol, s.M)

    # determine the number of link to be removed
    m = min(rand(1:u.var.m), u.var.max) # m-1 beads are removed
    jₘ = j₀ + m

    # compute the inital action
    w_initial = 0.0 # initial weight
    for j in j₀:(jₘ-1) # only these link contribute
        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        interaction_action!(w_initial, s, fpcycle, j)
    end

    # compute propagator
    Δk = prop_0(s.world[n].r[j₀, :], s.world[fpcycle(jₘ)].r[mod1(jₘ, s.M), :], m * s.τ, s.L, s.λ, s.dim)
    preamble = get_C(s.C) * u.var.m * s.N * s.M / Δk
    #metropolis question
    if metropolis(preamble * exp(-w_initial - s.μ * m * s.τ))
        acc = true
        s.Nz = s.N
        nw = ((jₘ - 1) ÷ s.M) + 1 # do we cross the boundary?

        # sidecases:
        # [1] Npol==1 && nw>1
        # [2] Npol==1 && nw==1
        # [3] Npol > 1 && nw == 1
        # [4] Npol ==2 && nw == 2
        # [5] Npol > 2 && nw ==2
        if Npol == 1 && nw > 1 # [1] Npol==1; nw>1
            r = Matrix{Union{Missing,Float64}}(missing, s.M, s.dim)
            V = Vector{Union{Missing,Float64}}(missing, s.M)
            bins = Vector{Union{Missing,Int64}}(missing, s.M)

            for j in mod1(jₘ, s.M):j₀
                r[j, :] = @view s.world[n].r[j, :]
                bins[j] = s.world[n].bins[j]
            end
            for j in mod1(jₘ, s.M):mod1(j₀ - 1, s.M)
                V[j] = s.world[n].V[j]
            end

            rm_nn!(s.nn, n, s.world[n].bins) # eerste removen is veiliger
            add_nn!(s.nn, length(s.world) + 1, bins)
            push!(s.world, Worm(r, V, bins, mod1(jₘ, s.M), j₀, 0))
            fixindexing!(s.world, s.nn, [n])
            deleteat!(s.world, n)
            s.N -= 1
            s.worms = 1
        else
            rt = Matrix{Union{Missing,Float64}}(missing, s.M, s.dim)
            Vt = Vector{Union{Missing,Float64}}(missing, s.M)
            bt = Vector{Union{Missing,Int64}}(missing, s.M)
            rh = Matrix{Union{Missing,Float64}}(missing, s.M, s.dim)
            Vh = Vector{Union{Missing,Float64}}(missing, s.M)
            bh = Vector{Union{Missing,Int64}}(missing, s.M)

            for j in 1:j₀
                rh[j, :] = s.world[n].r[j, :]
                bh[j] = s.world[n].bins[j]
            end
            # if j₀=1, Vh[j] all missing
            if j₀ != 1
                for j in 1:(j₀-1)
                    Vh[j] = s.world[n].V[j]
                end
            end

            for j in jₘ:(s.M*nw)
                rt[mod1(j, s.M), :] = s.world[fpcycle(j)].r[mod1(j, s.M), :]
                bt[mod1(j, s.M)] = s.world[fpcycle(j)].bins[mod1(j, s.M)]
                Vt[mod1(j, s.M)] = s.world[fpcycle(j)].V[mod1(j, s.M)]
            end

            s.worms = 2

            if Npol == 1 && nw == 1 # [2] Npol==1 && nw==1
                rm_nn!(s.nn, n, s.world[n].bins)
                # s.N + 2 because after fixindexing! will do -1
                add_nn!(s.nn, length(s.world) + 1, bt)
                push!(s.world, Worm(rt, Vt, bt, mod1(jₘ, s.M), 0, s.N + 2))
                add_nn!(s.nn, length(s.world) + 1, bh)
                push!(s.world, Worm(rh, Vh, bh, 0, j₀, 0))

                fixindexing!(s.world, s.nn, [n])
                deleteat!(s.world, n)
                s.N -= 1

            elseif Npol > 1 && nw == 1 # [3] Npol > 1 && nw == 1
                rm_nn!(s.nn, n, s.world[n].bins)

                add_nn!(s.nn, length(s.world) + 1, bt)
                push!(s.world, Worm(rt, Vt, bt, mod1(jₘ, s.M), 0, s.world[n].next))
                add_nn!(s.nn, length(s.world) + 1, bh)
                push!(s.world, Worm(rh, Vh, bh, 0, j₀, 0))
                # s.N+=2 because after fixindexing! will do -1
                s.world[cycle_findprev(s.world, n)].next = s.N + 2

                fixindexing!(s.world, s.nn, [n])
                deleteat!(s.world, n)
                s.N -= 1

            elseif Npol == 2 && nw == 2 # [4] Npol ==2 && nw == 2
                nnext = s.world[n].next

                rm_nn!(s.nn, n, s.world[n].bins)
                rm_nn!(s.nn, nnext, s.world[nnext].bins)
                # s.N+=2 because after fixindexing! will do -2
                add_nn!(s.nn, length(s.world) + 1, bt)
                push!(s.world, Worm(rt, Vt, bt, mod1(jₘ, s.M), 0, s.N + 2))
                add_nn!(s.nn, length(s.world) + 1, bh)
                push!(s.world, Worm(rh, Vh, bh, 0, j₀, 0))

                fixindexing!(s.world, s.nn, [n, nnext])
                deleteat!(s.world, sort([n, nnext]))
                s.N -= 2

            elseif Npol > 2 && nw == 2 # [5] Npol > 2 && nw ==2
                nnext = s.world[n].next
                nprev = cycle_findprev(s.world, n)

                rm_nn!(s.nn, n, s.world[n].bins)
                rm_nn!(s.nn, nnext, s.world[nnext].bins)

                add_nn!(s.nn, length(s.world) + 1, bt)
                push!(s.world, Worm(rt, Vt, bt, mod1(jₘ, s.M), 0, s.world[nnext].next))
                add_nn!(s.nn, length(s.world) + 1, bh)
                push!(s.world, Worm(rh, Vh, bh, 0, j₀, 0))

                fixindexing!(s.world, s.nn, [n, nnext])
                s.world[nprev].next = s.N
                deleteat!(s.world, sort([n, nnext]))
                s.N -= 2

            end
        end
    end
    queue!(u.counter_var, acc)
    return acc
end

"""
    CloseGC(s::System, slices::Int64, <keyword arguments>)

Initiates the Close struct where the variables for the Close update is stored. The Close update is one of the updates of the Worms Algorithm [1]. It closes particles if the amount of beads is smaller than a threshold(`slices`) when the simulation is in the G-sector and thereby transform the configuration the the Z-sector. It is completempary to the open update.

# Keyword arguments
- `minacc=0.4`: The lower bound accaptance used to change the number of slices that can be closed by the close-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices that can be closed by the close-function.

[1] 10.1103/PhysRevE.74.036701
"""

struct CloseGC <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function CloseGC(d::Threshold; adj=10, range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::CloseGC)(s::System)::Bool
    acc = false
    #only if config in G-sector
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end

    # find the tail and head beads
    # ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
    ntail = s.N + 1
    jₘ = s.world[ntail].tail # tail
    Npol, pol = subcycle(s.world, ntail)
    nhead = pol[Npol]
    j₀ = s.world[nhead].head  # head

    # (m-1) new bead will be sampled
    m = mod(jₘ - j₀, s.M) # mod not mod1!!!

    newp = jₘ - j₀ < 0 && nhead != ntail ? 2 : 1
    N′ = s.N + newp
    if !s.gc && N′ != s.Nz
        # error("The new particle number will be $(N′), however it must be $(s.Nz)")
        return acc
    end

    # if the number of beads to be closed is larger than the threshold --> reject
    # m = 0 because in the Open update one cannot create the inverse (detailed balance)
    if m > u.var.m || m == 0
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the new beads
    r′ = zeros(Float64, m + 1, s.dim)
    r′[begin, :] = s.world[nhead].r[j₀, :]
    r′[end, :] = s.world[ntail].r[jₘ, :]
    bool = hardspherelevy!(r′, s, j₀)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the new links and add them to the updated action `w_updated`
    w_updated = 0.0
    V′ = zeros(Float64, m) # new links
    for j in 1:m
        V′[j] = s.lnV(r′[j, :], r′[j+1, :])

        exceptions = j == 1 ? [nhead] : Int64[]
        interaction_action!(w_updated, s, r′, j, j₀, exceptions)
    end
    w_updated += sum(V′) # updated weight

    # compute the propagator
    Δk = prop_0(r′[1, :], r′[end, :], m * s.τ, s.L, s.λ, s.dim)
    preamble = Δk / (get_C(s.C) * u.var.m * N′ * s.M)
    # metropolis
    if metropolis(preamble * exp(w_updated + s.μ * m * s.τ))
        acc = true

        # amount of new particles
        newp = jₘ - j₀ < 0 && nhead != ntail ? 2 : 1

        if nhead == ntail
            rm_nn!(s.nn, nhead, s.world[nhead].bins)
        else
            rm_nn!(s.nn, nhead, s.world[nhead].bins)
            rm_nn!(s.nn, ntail, s.world[ntail].bins)
        end

        jj, ll = 1, 1
        for i in 1:newp
            # add beads and links to the world
            for j in jj:m
                rⱼ = @view r′[j, :]
                s.world[i == 1 ? nhead : ntail].r[mod1(j₀ + j - 1, s.M), :] = rⱼ
                s.world[i == 1 ? nhead : ntail].bins[mod1(j₀ + j - 1, s.M)] = bin(rⱼ, s.nbins, s.L)
                if mod1(j₀ + j - 1, s.M) == s.M && newp == 2
                    jj = j + 1
                    break
                end
            end
            for j in ll:m
                s.world[i == 1 ? nhead : ntail].V[mod1(j₀ + j - 1, s.M)] = V′[j]
                if mod1(j₀ + j - 1, s.M) == s.M && newp == 2
                    ll = j + 1
                    break
                end
            end

            # fill further up if jₘ - j₀ > 0 && newp==1
            if jₘ - j₀ > 0 && newp == 1
                @views for j in jₘ:s.M
                    s.world[nhead].r[j, :] = s.world[ntail].r[j, :]
                    s.world[nhead].bins[j] = s.world[ntail].bins[j]
                    s.world[nhead].V[j] = s.world[ntail].V[j]
                end
            end

            # sidecases
            if newp == 2
                if s.world[ntail].next == nhead
                    if i == 1
                        next = s.N + 1
                    elseif i == 2
                        next = s.N + 2
                    end
                else
                    if i == 1
                        next = s.N + 1
                        s.world[cycle_findprev(s.world, nhead)].next = s.N + 2
                    elseif i == 2
                        next = s.world[ntail].next
                    end
                end

            elseif newp == 1
                if iszero(s.world[ntail].next) || s.world[ntail].next == nhead
                    next = s.N + 1
                else
                    next = s.world[ntail].next
                    s.world[cycle_findprev(s.world, nhead)].next = s.N + 1
                end
            end

            add_nn!(s.nn, s.N + newp - i + 1, s.world[i == 1 ? nhead : ntail].bins)
            insert!(s.world, s.N + 1,
                Particle(
                    convert(Matrix{Float64}, s.world[i == 1 ? nhead : ntail].r),
                    convert(Vector{Float64}, s.world[i == 1 ? nhead : ntail].V),
                    s.world[i == 1 ? nhead : ntail].bins, next)
            )

            nhead += 1
            ntail += 1

            if jₘ - j₀ > 0 && i == 1
                deleteat!(s.world, nhead)
                deleteat!(s.world, ntail)
            else
                deleteat!(s.world, i == 1 ? nhead : ntail)

            end
        end
        s.worms = 0
        s.N += newp

    end

    queue!(u.counter_var, acc)
    return acc
end
