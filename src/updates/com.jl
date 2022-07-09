"""
    PolymerCenterOfMass(s::System, step::Int64, <keyword arguments>)

[1]
"""

struct PolymerCenterOfMass <: UpdateC
    counter::Counter
    var::Step
    counter_var::Counter

    function PolymerCenterOfMass(
        s::System,
        step::Float64;
        minstep = 1e-1,
        maxstep = s.L / 2,
        minacc = 0.4,
        maxacc = 0.6,
        adj = 10,
        range = 10_000
    )
        new(
            Counter(),
            Step(step, minstep, maxstep, minacc, maxacc),
            Counter(adj=adj, range=range)
        )
    end
end

# move the whole path of one polymer by a distance d
function (u::PolymerCenterOfMass)(s::System)::Bool
    # Bool to be returned which sasy if the update is accepted or not
    acc = false
    if s.N == 0
        return acc
    end

    Nlist = collect(1:s.N)
    if s.worms > 0
        push!(Nlist, s.N+1)
    end
    n′ = rand(Nlist)
    Npol, pol = subcycle(s.world, n′)

    # Compute the initial action w_initial; Since the relative positions are not changed their is no kinetic contribution.
    # Add contribution of potenial to every bead
    w_initial = 0.0

    # for every bead in the cycle
    for n in pol
        w_initial += sum(skipmissing(s.world[n].V))
         # Add contribution of the interaction between neigbouring beads
        for j in eachindex(skipmissing(s.world[n].r[:, 1]))
            interaction_action!(w_initial, s, _ -> n, j)
        end
    end

    # compute new positions of the beads r′ moved by a distance d wrt r
    # FIXME Vector of Matrix so that no disallowmissing
    r′ = Array{Union{Float64, Missing}}(missing, (s.M, s.dim, Npol))
    bool = move_polymer!(r′, s, u.var.size, pol)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # Now update the action
    w_updated = 0.0
    V′ = Array{Union{Float64, Missing}}(missing, (s.M, Npol))
    for (i, n) in pairs(pol)
        for j in eachindex(skipmissing(@view s.world[n].r[:, 1]))
            (nnext, mnext) = next(s.world, s.M, n, j)
            if iszero(nnext)
                continue
            end
            ii = nnext==n ? i : mod1(i+1, Npol)
            V′[j, i] = s.lnV(disallowmissing(r′[j, :, i]), disallowmissing(r′[mnext, :, ii]))

            interaction_action!(w_updated, s, r′, j, i, ii, mnext, Int64[n])
        end
    end
    w_updated += sum(skipmissing(V′))

    # Metropolis question
    if metropolis(exp(w_updated - w_initial))
        acc = true

        nlast = pol[Npol]
        mlast = isa(s.world[nlast], Worm) ?  s.world[nlast].head : 0

        for (i, n) in pairs(pol)
            for j in eachindex(skipmissing(@view s.world[n].r[:, 1]))
                if (n, j)!=(nlast, mlast)
                    s.world[n].V[j] = V′[j, i]
                end
                s.world[n].r[j, :] = view(r′, j, :, i)
                update_nn_bead!(s, n, j, disallowmissing(r′[j, :, i]))
            end
        end
    end

    queue!(u.counter_var, acc)
    return acc
end

"""
    SingleCenterOfMass(s::System, step::Int64, <keyword arguments>)

[1]
"""

struct SingleCenterOfMass <: UpdateC
    counter::Counter
    var::Step
    counter_var::Counter

    function SingleCenterOfMass(
        s::System,
        step::Float64;
        minstep = 1e-1,
        maxstep = s.L / 2,
        minacc = 0.4,
        maxacc = 0.6,
        adj = 10,
        range = 10_000
    )
        new(
            Counter(),
            Step(step, minstep, maxstep, minacc, maxacc),
            Counter(adj=adj, range=range)
        )
    end
end

# move the whole path of one polymer by a distance d
function (u::SingleCenterOfMass)(s::System)::Bool
    # Bool to be returned which sasy if the update is accepted or not
    acc = false
    if s.N == 0
        return acc
    end

    # It is inefficient to try to replace a polymer of size M*N > M
    singleparticles = Int64[]
    for i in 1:s.N
        if s.world[i].next == i
            push!(singleparticles, i)
        end
    end
    if s.worms > 0
        ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
        _, pol = subcycle(s.world, ntail)
        wlen = lengthpolymer(s.world, pol)
        if wlen <= s.M
            push!(singleparticles, ntail)
        end
    end
    if isempty(singleparticles)
        return acc
    end

    #chose a ramdom single worldline
    n1 = rand(singleparticles)
    Npol, pol = subcycle(s.world, n1)

    # Compute the initial action w_initial; Since the relative positions are not changed their is no kinetic contribution.
    # Add contribution of potenial to every bead
    w_initial = 0.0
    for n in pol    # for every bead in the cycle
        w_initial += sum(skipmissing(s.world[n].V))
        for j in eachindex(skipmissing(@view s.world[n].r[:, 1]))
            # Add contribution of the interaction between neigbouring beads
            interaction_action!(w_initial, s, _ -> n, j)
        end
    end

    # compute new positions of the beads r′ moved by a distance d wrt r
    # FIXME Vector of Matrix so that no disallowmissing
    r′ = Array{Union{Float64, Missing}}(missing, (s.M, s.dim, Npol))
    bool = move_polymer!(r′, s, u.var.size, pol)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # Now update the action
    w_updated = 0.0
    V′ = Array{Union{Float64, Missing}}(missing, (s.M, Npol))
    for (i, n) in pairs(pol)
        for j in eachindex(skipmissing(@view s.world[n].r[:, 1]))
            (nnext, mnext) = next(s.world, s.M, n, j)
            if iszero(nnext)
                continue
            end
            ii = nnext==n ? i : i+1
            V′[j, i] = s.lnV(disallowmissing(r′[j, :, i]), disallowmissing(r′[mnext, :, ii]))

            interaction_action!(w_updated, s, r′, j, i, ii, mnext, Int64[n])
        end
    end
    w_updated += sum(skipmissing(V′))

    # Metropolis question
    if metropolis(exp(w_updated - w_initial))
        acc = true

        nlast = pol[Npol]
        mlast = isa(s.world[nlast], Worm) ?  s.world[nlast].head : 0

        for (i, n) in pairs(pol)
            for j in eachindex(skipmissing(s.world[n].r[:, 1]))
                if (n, j)!=(nlast, mlast)
                    s.world[n].V[j] = V′[j, i]
                end
                s.world[n].r[j, :] = view(r′, j, :, i)
                update_nn_bead!(s, n, j, disallowmissing(r′[j, :, i]))
            end
        end
    end

    queue!(u.counter_var, acc)
    return acc
end
