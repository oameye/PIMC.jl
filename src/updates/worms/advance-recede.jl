"""
    Advance(s::System, slices::Int64, <keyword arguments>)

Initiates the Advance struct where the variables for the Advance update is stored. The Advance update is one of the updates of the Worms Algorithm [1]. It advances worms by adding a certain beads(`slices`). It is completempary to the Recede update.

# Keyword arguments
- `maxslices = s.N*s.M`: The maximum slices the advance-function can add to the worm. These variable can become arbitrary long. The default is set to the number of time slices of a worldline times the number of particle.
- `minacc=0.4`: The lower bound accaptance used to change the number of slices inserted by the advance-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices inserted by the advance-function.

[1] 10.1103/PhysRevE.74.036701
"""

struct Advance <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function Advance(
        d::Threshold;
        adj=10,
        range=100_000
    )
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::Advance)(s::System)::Bool
    acc = false
    # only if config in G-sector
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end

    # find tail and head beads
    ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
    Npol, pol = subcycle(s.world, ntail)
    nhead = pol[Npol]
    j₀ = s.world[nhead].head

    # determine slice to advance to
    m = rand(1:u.var.m)
    jₘ = j₀ + m

    # importance sampling in the canonical ensemble
    if !s.gc
        #head larger than tail
        hlt = s.world[nhead].head >= s.world[ntail].tail ? true : false
        if hlt
            if jₘ > s.M
                crossed = mod(jₘ, s.M) >= s.world[ntail].tail ? true : false
            else
                crossed = false
            end
        else
            if jₘ > s.M
                crossed = true
            else
                crossed = jₘ >= s.world[ntail].tail ? true : false
            end
        end
        if crossed
            return acc
        end
    end

    if ntail == nhead
        nprevhead = 0
    else
        nprevhead = cycle_findprev(s.world, nhead)
    end

    # generate new beads
    r′ = zeros(Float64, m + 1, s.dim)
    r′[begin, :] = s.world[nhead].r[j₀, :] # Masha
    r′[end, :] = gaussian_sample_hardsphere(s, r′[begin, :], jₘ, m)
    if all(r′[end, :] .== 0.0)
        queue!(u.counter_var, acc)
        return acc
    end

    bool = hardspherelevy!(r′, s, j₀)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # copmute the relevant links and the action
    w_updated = 0.0
    V′ = zeros(Float64, m) # new links
    for j in 1:m
        V′[j] = s.lnV(r′[j, :], r′[j+1, :])

        #filter nhead as otherwise we have it as nn neigbour at j_0
        exceptions = j == 1 ? [nhead] : Int64[]
        interaction_action!(w_updated, s, r′, j, j₀, exceptions)
    end
    # updated weight; there is no initial weight as only beads are added, not changed
    w_updated += sum(V′) # updated weight

    # metropolis question
    if metropolis(exp(w_updated + s.μ * m * s.τ))
        # check_bead_link(s)
        # check_nn(s)
        acc = true
        # times one goes over β-boundary
        nw = ((jₘ - 1) ÷ s.M) + 1

        # if this case we make two tangling ends
        if ntail == nhead && nw > 1
            s.worms = 2
        end

        rm_nn!(s.nn, nhead, s.world[nhead].bins)

        #ceate new worldlines
        jj, ll = 1, 1
        for i in 1:nw
            if i > 1
                s.world[nhead].r .= missing
                s.world[nhead].V .= missing
                s.world[nhead].bins .= missing
            end

            for j in jj:(m+1)
                rⱼ = @view r′[j, :]
                s.world[nhead].r[mod1(j₀ + j - 1, s.M), :] = rⱼ
                s.world[nhead].bins[mod1(j₀ + j - 1, s.M)] = bin(rⱼ, s.nbins, s.L)
                if mod1(j₀ + j - 1, s.M) == s.M
                    jj = j + 1
                    break
                end
            end
            for j in ll:m
                s.world[nhead].V[mod1(j₀ + j - 1, s.M)] = V′[j]
                if mod1(j₀ + j - 1, s.M) == s.M
                    ll = j + 1
                    break
                end
            end

            #if nw > 2 other particles get moved and must be renamed in nn

            # go through sidecases
            if nw == i
                s.world[nhead].head = mod1(jₘ, s.M)
                add_nn!(s.nn, nhead, s.world[nhead].bins)
            elseif i == nw - 1
                if i == 1 && !iszero(nprevhead)
                    s.world[nprevhead].next = s.N + nw - 1

                    add_nn!(s.nn, s.N + 1, s.world[nhead].bins)
                    for n in s.N+1:(length(s.world)-1)
                        fixindexing_nn!(s.nn, n, n + 1, s.world[n].bins)
                    end
                    # fixindexing_nn!(s.nn, ntail, ntail+1, s.world[ntail].bins)

                    insert!(s.world, s.N + 1, Particle(copy(s.world[nhead].r), copy(s.world[nhead].V), copy(s.world[nhead].bins), nhead + 1))
                    ntail += 1
                    nhead += 1
                elseif i == 1 && iszero(nprevhead)
                    add_nn!(s.nn, length(s.world), s.world[nhead].bins)

                    insert!(s.world, length(s.world), Worm(copy(s.world[nhead].r), copy(s.world[nhead].V), copy(s.world[nhead].bins), s.world[nhead].tail, 0, nhead + 1))

                    nhead += 1
                    s.world[nhead].tail = 0
                else
                    add_nn!(s.nn, s.N + 1, s.world[nhead].bins)
                    for n in s.N+1:(length(s.world)-1)
                        fixindexing_nn!(s.nn, n, n + 1, s.world[n].bins)
                    end

                    insert!(s.world, s.N + 1, Particle(copy(s.world[nhead].r), copy(s.world[nhead].V), copy(s.world[nhead].bins), nhead + 1))
                    ntail += 1
                    nhead += 1
                end
            else
                if i == 1 && !iszero(nprevhead)
                    s.world[nprevhead].next = s.N + nw - 1

                    add_nn!(s.nn, s.N + 1, s.world[nhead].bins)
                    for n in s.N+1:(length(s.world)-1)
                        fixindexing_nn!(s.nn, n, n + 1, s.world[n].bins)
                    end

                    insert!(s.world, s.N + 1, Particle(copy(s.world[nhead].r), copy(s.world[nhead].V), copy(s.world[nhead].bins), s.N + nw - i - 1))
                    ntail += 1
                    nhead += 1
                elseif i == 1 && iszero(nprevhead)
                    add_nn!(s.nn, length(s.world), s.world[nhead].bins)

                    insert!(s.world, length(s.world), Worm(copy(s.world[nhead].r), copy(s.world[nhead].V), copy(s.world[nhead].bins), s.world[nhead].tail, 0, s.N + nw - i - 1))

                    nhead += 1
                    s.world[nhead].tail = 0
                else
                    add_nn!(s.nn, s.N + 1, s.world[nhead].bins)
                    for n in s.N+1:(length(s.world)-1)
                        fixindexing_nn!(s.nn, n, n + 1, s.world[n].bins)
                    end

                    insert!(s.world, s.N + 1, Particle(copy(s.world[nhead].r), copy(s.world[nhead].V), copy(s.world[nhead].bins), s.N + nw - i - 1))
                    tail += 1
                    nhead += 1
                end
            end
        end
        s.N = iszero(nprevhead) ? max(s.N + nw - 2, s.N) : max(s.N + nw - 1, s.N)

    end
    queue!(u.counter_var, acc)
    return acc
end

"""
    Recede(s::System, slices::Int64, <keyword arguments>)

Initiates the Recede struct where the variables for the Recede update is stored. The Recede update is one of the updates of the Worms Algorithm [1]. It recedes worms by deleting a certain beads(`slices`). It is completempary to the Advance update.

# Keyword arguments
- `maxslices = s.N*s.M`: The maximum slices the Recede-function can delete to the worm. These variable can become arbitrary big. The default is set to the number of time slices of a worldline times the number of particle.
- `minacc=0.4`: The lower bound accaptance used to change the number of slices deleted by the Recede-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices deleted by the Recede-function.

[1] 10.1103/PhysRevE.74.036701
"""

struct Recede <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function Recede(
        d::Threshold;
        adj=10,
        range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::Recede)(s::System)::Bool
    acc = false
    # only in the G-sector
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end

    # find tail and head beads
    ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
    Npol, pol = subcycle(s.world, ntail)
    nhead = pol[Npol]
    jₘ = s.world[nhead].head

    # choose timeslice
    m = rand(1:u.var.m)
    # new head time slice
    j₀ = jₘ - m

    # importance sampling in the canonical ensemble
    if !s.gc
        #head larger than tail
        hlt = s.world[nhead].head >= s.world[ntail].tail ? true : false
        if hlt
            if j₀ > 0
                crossed = j₀ >= s.world[ntail].tail ? false : true
            else
                crossed = true
            end
        else
            if j₀ > 0
                crossed = false
            else
                crossed = mod1(j₀, s.M) >= s.world[ntail].tail ? false : true
            end
        end
        if crossed
            return acc
        end
    end

    # if more beads will be deleted than there are --> reject
    if m >= lengthpolymer(s.world, pol)
        queue!(u.counter_var, acc)
        return acc
    end

    # number of worldline to delete
    i = j₀ > 0 ? i = 0 : abs(((j₀) ÷ s.M) - 1)
    # worldlines to be removed
    pol = pol[end-i:end]

    fpcycle(j::Int64)::Int64 = pol[end+floor(Int64, (j - 1) / s.M)]

    # compute the inital action; there is no updated action because beads are removed
    w_initial = 0.0
    for j in (jₘ-m):jₘ-1
        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        interaction_action!(w_initial, s, fpcycle, j)
    end

    # metropolis question
    if metropolis(exp(-w_initial - s.μ * m * s.τ))
        acc = true

        # create new head
        rm_nn!(s.nn, nhead, s.world[nhead].bins)
        for j in 1:mod1(j₀, s.M)
            s.world[nhead].r[j, :] = s.world[pol[begin]].r[j, :]
            s.world[nhead].bins[j] = s.world[pol[begin]].bins[j]
        end
        s.world[nhead].r[(mod1(j₀, s.M)+1):end, :] .= missing
        s.world[nhead].bins[(mod1(j₀, s.M)+1):end] .= missing
        add_nn!(s.nn, nhead, s.world[nhead].bins)
        for j in 1:(mod1(j₀, s.M)-1)
            s.world[nhead].V[j] = s.world[pol[begin]].V[j]
        end
        s.world[nhead].V[mod1(j₀, s.M):end] .= missing
        s.world[nhead].head = mod1(j₀, s.M)

        # sidecases
        if length(pol) > 1 && isa(s.world[pol[begin]], Worm)
            s.world[nhead].tail = s.world[pol[begin]].tail
            s.worms = 1
        end
        if !isa(s.world[pol[begin]], Worm)
            s.world[cycle_findprev(s.world, pol[begin])].next = nhead
        end

        # delete worldlines
        if !isempty(pol[begin:end-1])
            rm_nn!(s.nn, pol[begin:end-1], s.world)
            fixindexing!(s.world, s.nn, pol[begin:end-1])
            s.N -= count(n -> !isa(s.world[n], Worm), pol[begin:end-1])
            deleteat!(s.world, sort(pol[begin:end-1]))
        end

    end
    queue!(u.counter_var, acc)
    return acc
end
