"""
    OpenTailC(s::System, slices::Int64, <keyword arguments>)

Initiates the Open struct where the variables for the Open update is stored. The Open update is one of the updates of the Worms Algorithm [1]. It opens particles with a certain amount of beads(`slices`) when the simulation is in the Z-sector and thereby transform the configuration the the G-sector. It is completempary to the Close update.

# Keyword arguments
- `maxslices = s.M-1`: The maximum slices the Open-function can Open into the simulation.
- `minacc=0.4`: The lower bound accaptance used to change the number of slices opened by the open-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices opened by the Open-function.

[1]
"""

struct OpenTailC <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function OpenTailC(d::Threshold; adj=10, range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::OpenTailC)(s::System)::Bool
    acc = false
    # only if config in Z-sector
    if s.worms > 0 || s.worm_algorithm == false
        return acc
    end

    # choose a bead random
    j₀, n = rand(1:s.M), rand(1:s.N)

    # determine the polymer
    Npol, pol = subcycle(s.world, n)
    fpcycle(j::Int64)::Int64 = pcycle(j, pol, Npol, s.M)

    # determine the number of link to be changed
    # FIXME
    # m = min(rand(1:u.var.m), u.var.max) # m beads are moved
    m = min(u.var.m, u.var.max)
    jₘ = j₀ + m

    # generate new beads
    r′ = zeros(Float64, m + 1, s.dim)
    r′[end, :] = s.world[fpcycle(jₘ)].r[mod1(jₘ, s.M), :]
    r′[begin, :] = gaussian_sample_hardsphere(s, r′[end, :], j₀, m) # no exceptions as (n, j₀) will remain for head
    if all(r′[begin, :] .== 0.0)
        queue!(u.counter_var, acc)
        return acc
    end
    bool = hardspherelevy!(r′, s, j₀, fpcycle)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the relevant links and the action
    w_initial = 0.0 # initial weight
    w_updated = 0.0
    V′ = zeros(Float64, m) # new links
    for j in j₀:(jₘ-1) # only these link contribute
        j′ =  j - j₀ + 1

        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        V′[j′] = s.lnV(r′[j′, :], r′[j′+1, :])

        interaction_action!(w_initial, s, fpcycle, j)
        interaction_action!(w_updated, s, r′, j′, j₀, [fpcycle(j)])
    end
    w_updated += sum(V′) # updated weight

    # compute propagator
    Δk = prop_0(r′[begin, :], r′[end, :], m * s.τ, s.L, s.λ, s.dim)
    # preamble = get_C(s.C) * s.N * s.M / Δk
    preamble = get_C(s.C) / Δk
    #metropolis question
    if metropolis(preamble * exp(w_updated - w_initial))
        acc = true
        s.Nz = s.N

        rt = Matrix{Union{Missing,Float64}}(missing, s.M, s.dim)
        Vt = Vector{Union{Missing,Float64}}(missing, s.M)
        bt = Vector{Union{Missing,Int64}}(missing, s.M)
        rh = Matrix{Union{Missing,Float64}}(missing, s.M, s.dim)
        Vh = Vector{Union{Missing,Float64}}(missing, s.M)
        bh = Vector{Union{Missing,Int64}}(missing, s.M)

        for j in 1:j₀
            rh[j, :] = s.world[n].r[j, :]
            bh[j] = s.world[n].bins[j]
            if j != j₀
                Vh[j] = s.world[n].V[j]
            end
        end
        for j in 1:m
            if j₀+j-1 <= s.M
            rt[mod1(j₀+j-1, s.M), :] = r′[j,:]
            bt[mod1(j₀+j-1, s.M)] = bin(r′[j,:], s.nbins, s.L)
            Vt[mod1(j₀+j-1, s.M)] = V′[j]
            elseif Npol > 1
                s.world[fpcycle(j₀+j-1)].r[mod1(j₀+j-1, s.M), :] = r′[j,:]
                # s.world[fpcycle(j₀+j-1)].bins[mod1(j₀+j-1, s.M)] = bin(r′[j,:], s.nbins, s.L)
                s.world[fpcycle(j₀+j-1)].V[mod1(j₀+j-1, s.M)] = V′[j]
                update_nn_bead!(s, fpcycle(j₀+j-1), mod1(j₀+j-1, s.M),  r′[j,:])
            else
                rh[mod1(j₀+j-1, s.M), :] = r′[j,:]
                bh[mod1(j₀+j-1, s.M)] = bin(r′[j,:], s.nbins, s.L)
                Vh[mod1(j₀+j-1, s.M)] =  V′[j]
            end
        end
        for j in jₘ:s.M
            rt[j, :] = s.world[n].r[j, :]
            bt[j] =  s.world[n].bins[j]
            Vt[j] = s.world[n].V[j]
        end

        # sidecases:
        # [1] Npol==1 && nw>1
        # [2] Npol>1 && nw==1
        tailnext = Npol > 1 ? s.world[n].next : s.N + 2 # s.N + 2 because after fixindexing! will do -1
        prevn = cycle_findprev(s.world, n)

        rm_nn!(s.nn, n, s.world[n].bins)
        add_nn!(s.nn, length(s.world) + 1, bt)
        push!(s.world, Worm(rt, Vt, bt, j₀, 0, tailnext))
        add_nn!(s.nn, length(s.world) + 1, bh)
        push!(s.world, Worm(rh, Vh, bh, 0, j₀, 0))

        s.world[prevn].next = length(s.world)

        fixindexing!(s.world, s.nn, [n])
        deleteat!(s.world, n)
        s.N -= 1
        s.worms = 2

    end

    queue!(u.counter_var, acc)
    return acc
end

"""
    OpenHeadC(s::System, slices::Int64, <keyword arguments>)

Initiates the Open struct where the variables for the Open update is stored. The Open update is one of the updates of the Worms Algorithm [1]. It opens particles with a certain amount of beads(`slices`) when the simulation is in the Z-sector and thereby transform the configuration the the G-sector. It is completempary to the Close update.

# Keyword arguments
- `maxslices = s.M-1`: The maximum slices the Open-function can Open into the simulation.
- `minacc=0.4`: The lower bound accaptance used to change the number of slices opened by the open-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices opened by the Open-function.

[1]
"""

struct OpenHeadC <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function OpenHeadC(d::Threshold; adj=10, range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::OpenHeadC)(s::System)::Bool
    acc = false
    # only if config in Z-sector
    if s.worms > 0 || s.worm_algorithm == false
        return acc
    end

    # choose a bead random
    jₘ, n = rand(1:s.M), rand(1:s.N)

    # determine the polymer
    Npol, pol = subcycle(s.world, n)
    fpcycle(j::Int64)::Int64 = pcycle(j, pol, Npol, s.M)

    # determine the number of link to be changed
    m = min(u.var.m, u.var.max)
    j₀ = jₘ - m

    # generate new beads
    r′ = zeros(Float64, m + 1, s.dim)
    r′[begin, :] = s.world[fpcycle(j₀)].r[mod1(j₀, s.M), :]
    r′[end, :] = gaussian_sample_hardsphere(s, r′[begin, :], jₘ, m) # no exceptions as (n, jₘ) will remain for head
    if all(r′[end, :] .== 0.0)
        queue!(u.counter_var, acc)
        return acc
    end
    bool = hardspherelevy!(r′, s, j₀, fpcycle)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the relevant links and the action
    w_initial = 0.0 # initial weight
    w_updated = 0.0
    V′ = zeros(Float64, m) # new links
    for j in j₀:(jₘ-1) # only these link contribute
        j′ = j - j₀ + 1

        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        V′[j′] = s.lnV(r′[j′, :], r′[j′+1, :])

        interaction_action!(w_initial, s, fpcycle, j)
        interaction_action!(w_updated, s, r′, j′, j₀, [fpcycle(j)])
    end
    w_updated += sum(V′) # updated weight

    # compute propagator
    Δk = prop_0(r′[begin, :], r′[end, :], m * s.τ, s.L, s.λ, s.dim)
    # preamble = get_C(s.C) * s.N * s.M / Δk
    preamble = get_C(s.C) / Δk
    #metropolis question
    if metropolis(preamble * exp(w_updated - w_initial))
        acc = true
        s.Nz = s.N

        rt = Matrix{Union{Missing,Float64}}(missing, s.M, s.dim)
        Vt = Vector{Union{Missing,Float64}}(missing, s.M)
        bt = Vector{Union{Missing,Int64}}(missing, s.M)
        rh = Matrix{Union{Missing,Float64}}(missing, s.M, s.dim)
        Vh = Vector{Union{Missing,Float64}}(missing, s.M)
        bh = Vector{Union{Missing,Int64}}(missing, s.M)

        for j in 1:j₀
            rh[j, :] = s.world[n].r[j, :]
            bh[j] = s.world[n].bins[j]
            if j != j₀
                Vh[j] = s.world[n].V[j]
            end
        end
        for j′ in 1:(m+1)
            j = j₀ + j′ - 1
            if j > 0
                rh[mod1(j, s.M), :] = r′[j′,:]
                bh[mod1(j, s.M)] = bin(r′[j′,:], s.nbins, s.L)
                if j != jₘ
                    Vh[mod1(j, s.M)] = V′[j′]
                end
            elseif Npol > 1
                s.world[fpcycle(j)].r[mod1(j, s.M), :] = r′[j′,:]
                # s.world[fpcycle(j₀+j-1)].bins[mod1(j₀+j-1, s.M)] = bin(r′[j,:], s.nbins, s.L)
                s.world[fpcycle(j)].V[mod1(j, s.M)] = V′[j′]
                update_nn_bead!(s, fpcycle(j), mod1(j, s.M),  r′[j′,:])
            else
                rt[mod1(j, s.M), :] = r′[j′,:]
                bt[mod1(j, s.M)] = bin(r′[j′,:], s.nbins, s.L)
                Vt[mod1(j, s.M)] =  V′[j′]
            end
        end
        jmax = Npol < 2 && j₀ < 0 ? mod1(j₀, s.M) : s.M
        for j in jₘ:jmax
            rt[j, :] = s.world[n].r[j, :]
            bt[j] =  s.world[n].bins[j]
            if ismissing(Vt[j])
                Vt[j] = s.world[n].V[j]
            end
        end

        # sidecases:
        # [1] Npol==1 && nw>1
        # [2] Npol>1 && nw==1
        tailnext = Npol > 1 ? s.world[n].next : s.N + 2 # s.N + 2 because after fixindexing! will do -1
        prevn = cycle_findprev(s.world, n)

        rm_nn!(s.nn, n, s.world[n].bins)
        add_nn!(s.nn, length(s.world) + 1, bt)
        push!(s.world, Worm(rt, Vt, bt, jₘ, 0, tailnext))
        add_nn!(s.nn, length(s.world) + 1, bh)
        push!(s.world, Worm(rh, Vh, bh, 0, jₘ, 0))

        s.world[prevn].next = length(s.world)

        fixindexing!(s.world, s.nn, [n])
        deleteat!(s.world, n)
        s.N -= 1
        s.worms = 2

    end

    queue!(u.counter_var, acc)
    return acc
end


"""
    CloseTailC(s::System, slices::Int64, <keyword arguments>)

Initiates the Close struct where the variables for the Close update is stored. The Close update is one of the updates of the Worms Algorithm [1]. It closes particles if the amount of beads is smaller than a threshold(`slices`) when the simulation is in the G-sector and thereby transform the configuration the the Z-sector. It is completempary to the open update.

# Keyword arguments
- `minacc=0.4`: The lower bound accaptance used to change the number of slices that can be closed by the close-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices that can be closed by the close-function.

[1]
"""
struct CloseTailC <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function CloseTailC(d::Threshold; adj=10, range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::CloseTailC)(s::System)::Bool
    acc = false
    # only if config in G-sector
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end

    # find the tail and head beads
    # ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
    ntail = s.N + 1
    jₘ = s.world[ntail].tail # tail
    Npol, pol = subcycle(s.world, ntail)
    fpcycle(j::Int64)::Int64 = pcycle(j, pol, Npol, s.M)
    nhead = pol[Npol]
    j₀ = s.world[nhead].head  # head

    if j₀ != jₘ
        error("In the canonical ensemble the timeslice of the head and worm  are not equal, i.e., j₀ != jₘ.")
    end

    # FIXME
    # m = min(rand(1:u.var.m), u.var.max) # m beads are moved
    m = min(u.var.m, u.var.max)
    jₘ = j₀ + m

    # generate new beads
    r′ = zeros(Float64, m + 1, s.dim)
    r′[begin, :] = s.world[nhead].r[j₀, :]
    r′[end, :] = s.world[fpcycle(jₘ)].r[mod1(jₘ, s.M), :]

    bool = hardspherelevy!(r′, s, j₀, fpcycle)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the relevant links and the action
    w_initial = 0.0 # initial weight
    w_updated = 0.0
    V′ = zeros(Float64, m) # new links
    for j in j₀:(jₘ-1) # only these link contribute
        j′ =  j - j₀ + 1

        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        interaction_action!(w_initial, s, fpcycle, j)

        V′[j′] = s.lnV(r′[j′, :], r′[j′+1, :])
        interaction_action!(w_updated, s, r′, j′, j₀, [fpcycle(j)])
    end
    w_updated += sum(V′) # updated weight

    # compute propagator
    Δk =  prop_0(r′[begin, :], r′[end, :], m * s.τ, s.L, s.λ, s.dim)
    # preamble = Δk / (get_C(s.C) * s.Nz * s.M)
    preamble = Δk / get_C(s.C)
    # metropolis
    if metropolis(preamble * exp(w_updated - w_initial))
        acc = true

        rm_nn!(s.nn, nhead, s.world[nhead].bins)
        rm_nn!(s.nn, ntail, s.world[ntail].bins)

        for j in 1:m
            if j₀ + j - 1 <= s.M
                s.world[nhead].r[j₀ + j - 1, :] =  r′[j, :]
                s.world[nhead].bins[j₀ + j - 1] = bin(r′[j, :], s.nbins, s.L)
                s.world[nhead].V[j₀ + j - 1] = V′[j]
            elseif Npol > 2
                s.world[fpcycle(j₀+j-1)].r[mod1(j₀+j-1, s.M), :] = r′[j,:]
                # s.world[fpcycle(j₀+j-1)].bins[mod1(j₀+j-1, s.M)] = bin(r′[j,:], s.nbins, s.L)
                s.world[fpcycle(j₀+j-1)].V[mod1(j₀+j-1, s.M)] = V′[j]
                update_nn_bead!(s, fpcycle(j₀+j-1), mod1(j₀+j-1, s.M),  r′[j,:])
            else
                s.world[nhead].r[mod1(j₀ + j - 1, s.M), :] =  r′[j, :]
                s.world[nhead].bins[mod1(j₀ + j - 1, s.M)] = bin(r′[j, :], s.nbins, s.L)
                s.world[nhead].V[mod1(j₀ + j - 1, s.M)] = V′[j]
            end
        end
        for j in jₘ:s.M
            s.world[nhead].r[j, :] = s.world[ntail].r[j, :]
            s.world[nhead].bins[j] =  s.world[ntail].bins[j]
            s.world[nhead].V[j] = s.world[ntail].V[j]
        end

        nexttail = Npol > 2 ? s.world[ntail].next : s.N + 1
        prevhead = Npol > 2 ? cycle_findprev(s.world, nhead) : ntail

        add_nn!(s.nn, s.N + 1, s.world[nhead].bins)
        s.world[prevhead].next = s.N + 1
        insert!(s.world, s.N + 1,
            Particle(
                convert(Matrix{Float64}, s.world[nhead].r),
                convert(Vector{Float64}, s.world[nhead].V),
                s.world[nhead].bins, nexttail
                )
        )

        deleteat!(s.world, nhead+1)
        deleteat!(s.world, ntail+1)

        s.worms = 0
        s.N += 1
    end

    queue!(u.counter_var, acc)
    return acc
end

"""
    CloseHeadC(s::System, slices::Int64, <keyword arguments>)

Initiates the Close struct where the variables for the Close update is stored. The Close update is one of the updates of the Worms Algorithm [1]. It closes particles if the amount of beads is smaller than a threshold(`slices`) when the simulation is in the G-sector and thereby transform the configuration the the Z-sector. It is completempary to the open update.

# Keyword arguments
- `minacc=0.4`: The lower bound accaptance used to change the number of slices that can be closed by the close-function.
- `maxacc=0.6`: The upper bound accaptance used to change the number of slices that can be closed by the close-function.

[1]
"""
struct CloseHeadC <: UpdateWorm
    counter::Counter
    var::Threshold
    counter_var::Counter

    function CloseHeadC(d::Threshold; adj=10, range=100_000)
        new(Counter(), d, Counter(adj=adj, range=range))
    end
end

function (u::CloseHeadC)(s::System)::Bool
    acc = false
    # only if config in G-sector
    if s.worms == 0 || s.worm_algorithm == false
        return acc
    end

    # find the tail and head beads
    # ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
    ntail = s.N + 1
    j₀ = s.world[ntail].tail # tail
    Npol, pol = subcycle(s.world, ntail)
    fpcycle(j::Int64)::Int64 = pcycle(j, pol, Npol, s.M)
    nhead = pol[Npol]
    jₘ = s.world[nhead].head  # head

    if j₀ != jₘ
        error("In the canonical ensemble the timeslice of the head and worm  are not equal, i.e., j₀ != jₘ.")
    end

    # FIXME
    # m = min(rand(1:u.var.m), u.var.max) # m beads are moved
    m = min(u.var.m, u.var.max)
    j₀ = jₘ - m
    pol = circshift(pol, 1)

    # generate new beads
    r′ = zeros(Float64, m + 1, s.dim)
    if ismissing(s.world[nhead].r[mod1(j₀, s.M), :])
        @info "debug"
    end
    r′[begin, :] = s.world[fpcycle(j₀)].r[mod1(j₀, s.M), :]
    r′[end, :] = s.world[ntail].r[mod1(jₘ, s.M), :]

    bool = hardspherelevy!(r′, s, j₀, fpcycle)
    if !bool
        queue!(u.counter_var, acc)
        return acc
    end

    # compute the relevant links and the action
    w_initial = 0.0 # initial weight
    w_updated = 0.0
    V′ = zeros(Float64, m) # new links
    for j in j₀:(jₘ-1) # only these link contribute
        j′ =  j - j₀ + 1

        w_initial += s.world[fpcycle(j)].V[mod1(j, s.M)]
        interaction_action!(w_initial, s, fpcycle, j)

        V′[j′] = s.lnV(r′[j′, :], r′[j′+1, :])
        interaction_action!(w_updated, s, r′, j′, j₀, [fpcycle(j)])
    end
    w_updated += sum(V′) # updated weight

    # compute propagator
    Δk =  prop_0(r′[begin, :], r′[end, :], m * s.τ, s.L, s.λ, s.dim)
    # preamble = Δk / (get_C(s.C) * s.Nz * s.M)
    preamble = Δk / get_C(s.C)
    # metropolis
    if metropolis(preamble * exp(w_updated - w_initial))
        acc = true

        rm_nn!(s.nn, nhead, s.world[nhead].bins)
        rm_nn!(s.nn, ntail, s.world[ntail].bins)

        jbool = false
        for j′ in 1:m
            j = j₀ + j′ - 1
            if j > 0
                s.world[nhead].r[mod1(j, s.M), :] =  r′[j′, :]
                s.world[nhead].bins[mod1(j, s.M)] = bin(r′[j′, :], s.nbins, s.L)
                s.world[nhead].V[mod1(j, s.M)] = V′[j′]
            elseif Npol > 2
                s.world[fpcycle(j)].r[mod1(j, s.M), :] = r′[j′,:]
                s.world[fpcycle(j)].V[mod1(j, s.M)] = V′[j′]
                update_nn_bead!(s, fpcycle(j), mod1(j, s.M),  r′[j′,:])
            else
                s.world[nhead].r[mod1(j, s.M), :] =  r′[j′, :]
                s.world[nhead].bins[mod1(j, s.M)] = bin(r′[j′, :], s.nbins, s.L)
                s.world[nhead].V[mod1(j, s.M)] = V′[j′]
                if j == s.M
                    jbool = true
                end
            end
        end
        jmax = Npol <= 2 && j₀ <= 0 ? mod1(j₀, s.M) : s.M
        for j in jₘ:jmax
            s.world[nhead].r[j, :] = s.world[ntail].r[j, :]
            s.world[nhead].bins[j] =  s.world[ntail].bins[j]
            if !jbool
                s.world[nhead].V[j] = s.world[ntail].V[j]
            end
        end

        nexttail = Npol > 2 ? s.world[ntail].next : s.N + 1
        prevhead = Npol > 2 ? cycle_findprev(s.world, nhead) : ntail

        add_nn!(s.nn, s.N + 1, s.world[nhead].bins)
        s.world[prevhead].next = s.N + 1
        insert!(s.world, s.N + 1,
            Particle(
                convert(Matrix{Float64}, s.world[nhead].r),
                convert(Vector{Float64}, s.world[nhead].V),
                s.world[nhead].bins, nexttail
                )
        )
        deleteat!(s.world, nhead+1)
        deleteat!(s.world, ntail+1)

        s.worms = 0
        s.N += 1
    end

    queue!(u.counter_var, acc)
    return acc
end