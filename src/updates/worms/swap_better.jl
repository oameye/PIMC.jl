"""
    Swap(s::System, slices::Int64, <keyword arguments>)


Initiates the Swap struct where the variables for the Swap update is stored.
The Swap update is a central move in the Worms algorithm [1] which proposes new configuration when the simulation is in
the g-sector. It swaps worms with particle worldlines and thereby enhanced the ability the propose configuration with a
high Winding number. The ladder is critical for simulating phenonema such as superfluidity.


# Arguments
- `minslices = 1`: The minimum slices the swapfunction can modify. Must be larger or equel one.
- `maxslices = s.M-1`: The maximum slices the swapfunction can modify. Must be smaller or equal the number of slices of a worldline minus one.
- `minacc=0.4`: The lower bound accaptance used to change the number of slices modified by the Swap function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices modified by the Swap function.

[1] 10.1103/PhysRevE.74.036701
"""
struct SwapC <: UpdateWorm
    counter::Counter
    var::NumbOfSlices
    counter_var::Counter

    function SwapC(d::NumbOfSlices; adj=10, range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::SwapC)(s::System)::Bool
    acc = false
    # Cannot use the Swap update when in the z-sector.
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end
    if s.N == 0
        return acc
    end

    # find the head and tail of the worms
    ntail = s.N + 1
    NpolW, polW = subcycle(s.world, ntail)
    nhead = polW[NpolW]

    # take the lower and upper bound for the to-be modified time slices
    if !isa(s.world[nhead], Worm)
        println(debugpermlist(s))
    end
    j₀ = s.world[nhead].head
    m = u.var.m # not rand!
    jₘ = mod1(j₀ + m, s.M)

    # sample a particle in the nearby bins to nhead at time slice jₘ
    # norm1 is needed for the accaptance probability in metropolis
    # only closed worldlines will be sampled
    n1, norm1 = sampleparticle_swap(s, m, nhead, j₀)

    # TODO a sidecase to implement later
    n1prev = cycle_findprev(s.world, n1)
    if nhead==ntail || n1prev == ntail
        return acc
    end

    # store the time slices which will be deleted in a list
    list = Vector{Tuple{Int64,Int64}}(undef, m + 1)
    for j in eachindex(list)
        list[j] = j == 1 ? (n1, jₘ) : prev(s.world, s.M, list[j-1]...)
    end

    # if the tail must be deleted to perform the update, reject the move
    for l in list
        # (ntail, s.world[ntail].tail) could be last item of the list
        if l[1] == 0 || (ntail, s.world[ntail].tail) == l
            queue!(u.counter_var, acc)
            return acc
        end
    end

    t2 = createKlist(s, m, list[m+1]...)
    norm2 = sum(t2)

    # compute the initial action
    w_initial = 0.0 # initial weight
    for j in 2:m+1 # only these link contribute
        w_initial += s.world[list[j][1]].V[list[j][2]]
        interaction_action!(w_initial, s, _ -> list[j][1], list[j][2])
    end

    # compute the new position with the levy construction
    r′ = zeros(Float64, m + 1, s.dim)
    r′[begin, :] = s.world[nhead].r[j₀, :]
    r′[end, :] = s.world[n1].r[jₘ, :]
    bool = hardspherelevy!(r′, s, j₀)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the updated action
    w_updated = 0.0
    V′ = zeros(Float64, m) # new links
    for j in eachindex(V′)
        V′[j] = s.lnV(r′[j, :], r′[j+1, :])

        exceptions=[list[j][1], nhead, n1] # n1 == list[j][1]?
        interaction_action!(w_updated, s, r′, j, j₀, exceptions)
    end
    w_updated += sum(V′) # updated weight

    # metroplis-hastings
    if metropolis((norm1 / norm2)*exp(w_updated - w_initial))
        acc = true
        pre_permlist = debugpermlist(s)

        # determine if the upper β-periodic boundary is crossed; yes nw=2 no nw=1
        nw = ((j₀ + m - 1) ÷ s.M) + 1
        # bcs n1 closed worldline prev/next never 0
        n1prev = cycle_findprev(s.world, n1)
        n1prevprev = cycle_findprev(s.world, n1prev)
        n1next = s.world[n1].next
        # can be 0, than we have nhead=ntail
        nheadprev = cycle_findprev(s.world, nhead)
        # if 0 ∈ [n1prev, n1prevprev, n1next, nheadprev]
        #     println("debug")
        # end

        rm_nn!(s.nn, nhead, s.world[nhead].bins)
        rm_nn!(s.nn, n1, s.world[n1].bins)
        rm_nn!(s.nn, n1prev, s.world[n1prev].bins)

        # new nhead
        newhead = nw==2 ? n1prev : n1
        for j = 1:j₀
            # if ismissing(s.world[nhead].r[j, :]) || ismissing(s.world[newhead].r[j, :])
            #     println("debug")
            # end
            s.world[nhead].r[j, :], s.world[newhead].r[j, :] = s.world[newhead].r[j, :], s.world[nhead].r[j, :]
            s.world[nhead].bins[j], s.world[newhead].bins[j] = s.world[newhead].bins[j], s.world[nhead].bins[j]
            if j !=j₀
                s.world[nhead].V[j], s.world[newhead].V[j] = s.world[newhead].V[j], s.world[nhead].V[j]
            end
        end
        # ^ beads newhead done
        # new particle worldline beads
        plist = nw==2 ? [n1prev, n1] : [n1]
        fpcycle(j) = plist[((j - 1) ÷ s.M) + 1]
        for j in j₀:(j₀+m-1) # jₘ is mod1ed
            j′ = j - j₀ + 1
            s.world[fpcycle(j)].r[mod1(j, s.M), :] = r′[j′, :]
            s.world[fpcycle(j)].bins[mod1(j, s.M)] = bin(r′[j′, :], s.nbins, s.L)
            s.world[fpcycle(j)].V[mod1(j, s.M)] = V′[j′]
        end

        add_nn!(s.nn, nhead, s.world[nhead].bins)
        add_nn!(s.nn, n1, s.world[n1].bins)
        add_nn!(s.nn, n1prev, s.world[n1prev].bins)

        # [3] nw == 1; ntail != nhead; n1 != ntail; n1 != s.world[n1].next
        # [4] nw == 1; ntail != nhead; n1 == s.world[n1].next
        # [7] nw == 2; ntail != nhead; n1 == n1prev
        # [8] nw == 2; ntail != nhead; n1 != n1prev

        # fix permutations
        if n1 ∉ polW
            if n1 == n1next
                s.world[n1].next = nhead
                s.world[nheadprev].next = n1
            else
                if nw == 1
                    s.world[n1].next = n1next
                    s.world[n1prev].next = nhead
                    s.world[nheadprev].next = n1
                else # nw == 2
                    s.world[nheadprev].next = n1prev
                    s.world[n1prev].next = n1
                    s.world[n1prevprev].next = nhead
                end
            end
        else
            if nw == 1
                if n1next == nhead
                    s.world[n1].next = n1
                    s.world[n1prev].next = nhead
                else
                    s.world[n1].next = n1next
                    s.world[n1prev].next = nhead
                    s.world[nheadprev].next = n1
                end
            else # nw == 2
                if n1next == nhead
                    s.world[n1].next = n1prev
                    s.world[n1prev].next = n1
                    s.world[n1prevprev].next = nhead
                else
                    s.world[n1].next = n1next
                    s.world[n1prev].next = n1
                    s.world[nheadprev].next = n1prev
                    s.world[n1prevprev].next = nhead
                end
            end
        end
        ntail = s.N + 1
        bool = subcycle(s, ntail, pre_permlist)
        if !bool
            error()
        end

        ntail = s.N + 1
        NpolW, polW = subcycle(s.world, ntail)
        nhead = polW[NpolW]
        if !isa(s.world[nhead], Worm)
            println(pre_permlist)
            println(debugpermlist(s))
        end
    end

    queue!(u.counter_var, acc)
    return acc
end
