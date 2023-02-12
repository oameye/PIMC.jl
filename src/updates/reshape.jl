"""
    ReshapeLinear(s::System, slices::Int64, <keyword arguments>)

[1]
"""

struct ReshapeLinear <: UpdateC
    counter::Counter
    var::NumbOfSlices
    counter_var::Counter

    function ReshapeLinear(
        s::System,
        slices::Int64;
        minslices = 2,
        maxslices = s.M - 2,
        minacc = 0.6,
        maxacc = 0.8,
        adj = 10,
        range = 10_000
    )

        new(
            Counter(),
            NumbOfSlices(min(slices, maxslices), minslices, maxslices, minacc, maxacc),
            Counter(adj = adj, range = range)
        )
    end
end

function (u::ReshapeLinear)(s::System)::Bool
    acc = false
    # Choose initial bead and determine its polymer
    if s.N == 0 && s.worms == 0
        return acc
    end
    n = rand(eachindex(s.world))
    j₀ = rand(collect(eachindex(skipmissing(s.world[n].r[:,1]))))
    Ncycle, cycle = subcycle(s.world, n)

    # make curry function wich gives particle index form time slice
    fpcycle(j::Int64) = pcycle(j, cycle, Ncycle, s.M)

    # determine
    m = min(u.var.maxslices, rand(2:u.var.m)) # smaller than s.M
    jₘ = j₀ + m

    # if s.worms > 0 && length(s.world) ∈ cycle[1:Ncycle]
    #     jhead = s.world[end].head
    #     mj = (s.M - j₀) + (Ncycle-2)*s.M + jhead
    #     if mj < m
    #         return acc
    #     end
    # end

    r′ = zeros(Float64, m + 1, s.dim)
    r′[1, :] = s.world[n].r[j₀, :]
    r′[m+1, :] = s.world[fpcycle(jₘ)].r[mod1(jₘ, s.M), :]

    bool = hardspherelevy!(r′, s, j₀, fpcycle)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    V′ = zeros(Float64, m)
    w_initial, w_updated = 0.0, 0.0 # initial weight andweight after update
    for j in j₀:(jₘ-1) # only these link contribute
        j′ = j - j₀ + 1

        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        interaction_action!(w_initial, s, fpcycle, j)

        V′[j′] = s.lnV(r′[j′, :], r′[j′+1, :])
        interaction_action!(w_updated, s, r′, j′, j₀, [fpcycle(j)])
    end
    w_updated += sum(V′)

    # Metropolis question
    if metropolis(exp(w_updated - w_initial))
        acc = true
        for j in eachindex(V′)
            s.world[fpcycle(j₀ + j - 1)].r[mod1(j₀ + j - 1, s.M), :] = @view r′[j, :]
            s.world[fpcycle(j₀ + j - 1)].V[mod1(j₀ + j - 1, s.M)] = V′[j]
            update_nn_bead!(s, fpcycle(j₀ + j - 1), mod1(j₀ + j - 1, s.M), view(r′, j, :))
        end
    end

    queue!(u.counter_var, acc)
    return acc
end

"""
    ReshapeSwapLinear(s::System, slices::Int64, <keyword arguments>)

[1]
"""

struct ReshapeSwapLinear <: UpdateC
    counter::Counter
    var::NumbOfSlices
    counter_var::Counter

    function ReshapeSwapLinear(
        s::System,
        slices::Int64;
        minslices = 2,
        maxslices = s.M - 2,
        minacc = 0.6,
        maxacc = 0.8,
        adj = 10,
        range = 10_000
    )

        new(
            Counter(),
            NumbOfSlices(min(slices, maxslices), minslices, maxslices, minacc, maxacc),
            Counter(adj = adj, range = range)
        )
    end
end

function (u::ReshapeSwapLinear)(s::System)::Bool
    acc = false
    if s.N <= 1
        return acc
    end # nothing to permute

    j₀ = rand(1:s.M)
    m = min(u.var.maxslices, rand(2:u.var.m))
    jₘ = j₀ + m

    n1, n2 = sampleparticles(s, j₀, m) # Choose two particles at time slice j₀
    if n1 == n2
        return acc
    end

    Npol1, pol1 = subcycle(s.world, n1)
    Npol2, pol2 = subcycle(s.world, n2)
    fpcycle1(j::Int64) = pcycle(j, pol1, Npol1, s.M)
    fpcycle2(j::Int64) = pcycle(j, pol2, Npol2, s.M)

    # if s.world[fpcycle2(jₘ)] isa Worm || s.world[fpcycle1(jₘ)] isa Worm
    #     return acc
    # end

    r1 = zeros(Float64, m + 1, s.dim)
    r2 = zeros(Float64, m + 1, s.dim)
    r1[1, :] = s.world[n1].r[j₀, :]
    r2[1, :] = s.world[n2].r[j₀, :]
    r1[m+1, :] = s.world[fpcycle2(jₘ)].r[mod1(jₘ, s.M), :]
    r2[m+1, :] = s.world[fpcycle1(jₘ)].r[mod1(jₘ, s.M), :]

    #pol and Npol for old positions in expandbox
    bool1 = hardspherelevy!(r1, s, j₀, fpcycle2)
    bool2 = hardspherelevy!(r2, s, j₀, fpcycle1)
    if !bool1 || !bool2
        queue!(u.counter_var, acc)
        return acc
    end

    w_initial = 0.0 # initial weight
    for j in j₀:(jₘ-1) # only these link contribute
        w_initial += s.world[fpcycle1(j)].V[mod1(j, s.M)] + s.world[fpcycle2(j)].V[mod1(j, s.M)]

        nnlist1 = find_nns(s, fpcycle1(j), mod1(j, s.M), exceptions = [fpcycle1(j)])
        nnlist2 = find_nns(s, fpcycle2(j), mod1(j, s.M), exceptions = [fpcycle2(j)])
        for nn in nnlist1
            nn_nextbead = next(s.world, s.M, nn, mod1(j, s.M))
            j_nextbead = next(s.world, s.M, fpcycle1(j), mod1(j, s.M))
            if iszero(nn_nextbead[1])
                continue
            end
            r = distance.(
                s.world[nn].r[mod1(j, s.M), :],
                s.world[fpcycle1(j)].r[mod1(j, s.M), :],
                s.L)
            r′ = distance.(
                s.world[nn_nextbead[1]].r[nn_nextbead[2], :],
                s.world[j_nextbead[1]].r[j_nextbead[2], :],
                s.L)
            w_initial += s.lnU(r, r′)
        end
        for nn in nnlist2
            nn_nextbead = next(s.world, s.M, nn, mod1(j, s.M))
            j_nextbead = next(s.world, s.M, fpcycle2(j), mod1(j, s.M))
            if iszero(nn_nextbead[1])
                continue
            end
            r = distance.(
                s.world[nn].r[mod1(j, s.M), :],
                s.world[fpcycle2(j)].r[mod1(j, s.M), :],
                s.L)
            r′ = distance.(
                s.world[nn_nextbead[1]].r[nn_nextbead[2], :],
                s.world[j_nextbead[1]].r[j_nextbead[2], :],
                s.L)
            w_initial += s.lnU(r, r′)
        end

    end

    V1, V2 = zeros(Float64, m), zeros(Float64, m) # new links
    for j in j₀:(jₘ-1)
        j′ = j - j₀  + 1
        V1[j′], V2[j′] = s.lnV(r1[j′, :], r1[j′+1, :]), s.lnV(r2[j′, :], r2[j′+1, :])

        # fpcycle1(j), fpcycle2(j) are in exceptions as they will not be in the update configuration
        nnlist1 = find_nns(s, r1[j′, :], mod1(j, s.M), exceptions = [fpcycle1(j), fpcycle2(j)])
        nnlist2 = find_nns(s, r2[j′, :], mod1(j, s.M), exceptions = [fpcycle1(j), fpcycle2(j)])
        for nn in nnlist1
            nn_nextbead = next(s.world, s.M, nn, mod1(j, s.M))
            if iszero(nn_nextbead[1])
                continue
            end
            r = distance.(
                s.world[nn].r[mod1(j, s.M), :],
                r1[j′, :],
                s.L)
            r′ = distance.(
                s.world[nn_nextbead[1]].r[nn_nextbead[2], :],
                r1[j′+1, :],
                s.L)
            w_initial += s.lnU(r, r′)
        end
        for nn in nnlist2
            nn_nextbead = next(s.world, s.M, nn, mod1(j, s.M))
            if iszero(nn_nextbead[1])
                continue
            end
            r = distance.(
                s.world[nn].r[mod1(j, s.M), :],
                r2[j′, :],
                s.L)
            r′ = distance.(
                s.world[nn_nextbead[1]].r[nn_nextbead[2], :],
                r2[j′+1, :],
                s.L)
            w_initial += s.lnU(r, r′)
        end
        # TODO add interactions between r1 and r2
    end
    w_updated = sum(V1)+ sum(V2) # updated weight


    # Metropolis question
    if metropolis(exp(w_updated - w_initial))
        acc =true

        rm_nn!(s.nn, pol1, s.world)
        rm_nn!(s.nn, pol2, s.world)

        #swap
        s.world[n1].next, s.world[n2].next = s.world[n2].next, s.world[n1].next

        Npol1, pol1 = subcycle(s.world, n1)
        Npol2, pol2 = subcycle(s.world, n2)

        for j in 2:(m+1)
            s.world[fpcycle1(j₀ + j - 1)].r[mod1(j₀ + j - 1, s.M), :] = r1[j, :]
            s.world[fpcycle1(j₀ + j - 1)].bins[mod1(j₀ + j - 1, s.M)] = bin(r1[j, :], s.nbins, s.L)
            s.world[fpcycle2(j₀ + j - 1)].r[mod1(j₀ + j - 1, s.M), :] = r2[j, :]
            s.world[fpcycle2(j₀ + j - 1)].bins[mod1(j₀ + j - 1, s.M)] = bin(r2[j, :], s.nbins, s.L)
        end
        for j in 1:m
            s.world[fpcycle1(j₀ + j - 1)].V[mod1(j₀ + j - 1, s.M)] = V1[j]
            s.world[fpcycle2(j₀ + j - 1)].V[mod1(j₀ + j - 1, s.M)] = V2[j]
        end
        if jₘ < s.M
            for j in (jₘ+1):s.M
                s.world[n1].r[j, :], s.world[n2].r[j, :] = s.world[n2].r[j, :], s.world[n1].r[j, :]
                s.world[n1].bins[j], s.world[n2].bins[j] = s.world[n2].bins[j], s.world[n1].bins[j]
                s.world[n1].V[j], s.world[n2].V[j] = s.world[n2].V[j], s.world[n1].V[j]
            end
        end

        add_nn!(s.nn, pol1, s.world)
        add_nn!(s.nn, pol2, s.world)
    end

    queue!(u.counter_var, acc)
    return acc
end
