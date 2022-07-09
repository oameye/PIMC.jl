struct MoveHead <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function MoveHead(
        d::Threshold;
        adj = 10,
        range = 100_000
    )
        new(Counter(), d, Counter(adj = adj, range = range))
    end
end

function (u::MoveHead)(s::System)::Bool
    acc = false
    # only if config in G-sector
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end

    # find tail and head beads
    # ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
    ntail = s.N + 1
    Npol, pol = subcycle(s.world, ntail)
    nhead = pol[Npol]
    jₘ = s.world[nhead].head

    # determine bead(= links) to be reshaped
    # m = rand(1:u.var.m)
    m = min(u.var.m, u.var.max)
    # if more beads will be reshaped than there are --> reject
    if m >= lengthpolymer(s.world, pol)
        queue!(u.counter_var, acc)
        return acc
    end
    # new head time slice
    j₀ = jₘ - m
    # number of worldlines involved
    i = j₀ > 0 ? i = 0 : abs((j₀ ÷ s.M) - 1)

    # worldlines to be reshaped
    pol = pol[end-i:end]

    # generate new beads
    r′ = zeros(Float64, m + 1, s.dim) # m+1 bead involved
    r′[begin, :] = s.world[pol[begin]].r[mod1(j₀, s.M), :] # Masha
    r′[end, :] = gaussian_sample_hardsphere(s, r′[begin, :], jₘ, m, exceptions = Int64[nhead])
    if all(r′[end, :] .== 0.0)
        queue!(u.counter_var, acc)
        return acc
    end

    fpcycle(j::Int64)::Int64 = pol[end+floor(Int64, (j - 1) / s.M)]
    bool = hardspherelevy!(r′, s, j₀, fpcycle)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the relevant links and the action
    w_updated = 0.0
    w_initial = 0.0 # initial weight
    V′ = zeros(Float64, m) # new links
    for j in j₀:(jₘ-1) # only these link contribute
        j′ = j - j₀ + 1
        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        interaction_action!(w_initial, s, fpcycle, j)

        V′[j′] = s.lnV(r′[j′, :], r′[j′+1, :])
        interaction_action!(w_updated, s, r′, j′, j₀, Int64[fpcycle(j)])
    end
    w_updated += sum(V′) # updated weight

    # metropolis question
    if metropolis(exp(w_updated - w_initial))
        acc = true
        for j in j₀:jₘ
            j′ = j - j₀ + 1
            s.world[fpcycle(j)].r[mod1(j, s.M), :] = r′[j′, :]
            update_nn_bead!(s, fpcycle(j), mod1(j, s.M), r′[j′, :])
            if j != jₘ
                s.world[fpcycle(j)].V[mod1(j, s.M)] = V′[j′]
            end
        end
    end

    queue!(u.counter_var, acc)
    return acc
end

struct MoveTail <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function MoveTail(
        d::Threshold;
        adj = 10,
        range = 100_000
    )
        new(Counter(), d, Counter(adj = adj, range = range))
    end
end

function (u::MoveTail)(s::System)::Bool
    acc = false
    # only if config in G-sector
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end

    # find tail and head beads
    # ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
    ntail = s.N + 1
    Npol, pol = subcycle(s.world, ntail)
    # nhead = pol[Npol]
    j₀ = s.world[ntail].tail

    # determine beads()links)= to reshape to
    # m = rand(1:u.var.m)
    m = min(u.var.m, u.var.max)

    # if more beads will be deleted than there are --> reject
    if m >= lengthpolymer(s.world, pol)
        queue!(u.counter_var, acc)
        return acc
    end
    # new head time slice
    jₘ = j₀ + m
    # number of worldlines involved
    i = 1 + (jₘ-1) ÷ s.M

    # worldlines to be reshaped
    pol = pol[begin:end-(Npol-i)]

    # generate new beads
    r′ = zeros(Float64, m + 1, s.dim)
    r′[end, :] = s.world[pol[end]].r[mod1(jₘ, s.M), :] # Masha
    r′[begin, :] = gaussian_sample_hardsphere(s, r′[end, :], j₀, m, exceptions=Int64[ntail])
    if all(r′[begin, :] .== 0.0)
        queue!(u.counter_var, acc)
        return acc
    end

    fpcycle(j::Int64)::Int64 = pol[begin+((j-1) ÷ s.M)]
    bool = hardspherelevy!(r′, s, j₀, fpcycle)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the relevant links and the action
    w_updated = 0.0
    w_initial = 0.0 # initial weight
    V′ = zeros(Float64, m) # new links
    for j in j₀:(jₘ-1) # only these link contribute
        j′ = j - j₀ + 1
        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        V′[j′] = s.lnV(r′[j′, :], r′[j′+1, :])

        interaction_action!(w_initial, s, fpcycle, j)
        interaction_action!(w_updated, s, r′, j′, j₀, Int64[fpcycle(j)])
    end
    w_updated += sum(V′) # updated weight

    # metropolis question
    if metropolis(exp(w_updated - w_initial))
        acc = true
        for j in j₀:jₘ
            j′ = j - j₀ + 1
            s.world[fpcycle(j)].r[mod1(j, s.M), :] = @view r′[j′, :]
            update_nn_bead!(s, fpcycle(j), mod1(j, s.M), @view r′[j′, :])
            if j != jₘ
                s.world[fpcycle(j)].V[mod1(j, s.M)] = V′[j′]
            end
        end
    end
    queue!(u.counter_var, acc)
    return acc
end