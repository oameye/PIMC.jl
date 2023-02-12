using Random, StatsBase

function metropolis(δ::Float64)
    (δ >= 1.0) || (δ > rand())
end
Base.@kwdef mutable struct Counter
    queue::Queue{Bool} = Queue{Bool}()
    tries::Int64 = 0
    range::Int64 = 100_000
    adj::Int64 = 10
end

# Some updates require a mutable step size; it is adjusted DURING an update
mutable struct Step
    size::Float64
    minstep::Float64
    maxstep::Float64
    minacc::Float64
    maxacc::Float64
end

function adjust!(step::Step, acceptance::Float64)
    if acceptance < step.minacc
        step.size *= 0.9
    elseif acceptance > step.maxacc
        step.size *= 1.1
    end

    # Keep within bounds
    step.size = max(step.minstep, step.size)
    step.size = min(step.maxstep, step.size)
end

mutable struct NumbOfSlices
    m::Int64
    minslices::Int64
    maxslices::Int64
    minacc::Float64
    maxacc::Float64
end

function adjust!(m::NumbOfSlices, acc::Float64)
    if acc < m.minacc
        m.m -= 1
    elseif acc > m.maxacc
        m.m += 1
    end

    # Keep within bounds
    m.m = max(m.minslices, m.m)
    m.m = min(m.maxslices, m.m)
end


function cycle_findprev(list::Vector{Worldline}, n::Int64)::Int64
    for (i, p) in pairs(list)
        if !iszero(p.next) && n == p.next
            return i
        end
    end
    return 0
end

function subcycle(p::Vector{Worldline}, n::Int64)::Tuple{Int64, Vector{Int64}}
    cycle = Int64[]
    push!(cycle, n)
    Ncycle = 1
    if !iszero(p[n].next) && p[n].next != n
        i = n
        ctr = 0
        while true
            ctr += 1
            i = p[i].next
            if iszero(i) || i == n
                break
            end
            Ncycle += 1
            push!(cycle, i)
            if ctr > length(p)
                error()
            end
        end
    end
    return Ncycle, cycle
end

function subcycle(s::System, n::Int64, preperm)
    p= s.world
    cycle = Int64[]
    push!(cycle, n)
    Ncycle = 1
    if !iszero(p[n].next) && p[n].next != n
        i = n
        ctr = 0
        while true
            ctr += 1
            i = p[i].next
            if iszero(i) || i == n
                break
            end
            Ncycle += 1
            push!(cycle, i)
            if ctr > length(p)*2
               println(preperm)
               println(debugpermlist(s))
               return false
            end
        end
    end
    return true
end

function pcycle(j::Int64, pol::Vector{Int64}, Npol::Int64, M::Int64)::Int64
    return pol[mod1(1 + floor(Int64, (j - 1) / M), Npol)]
end

# Builds a brownian bridge for a worldline in a time window
function levy!(r′::Matrix{Float64}, τ::Float64, L::Float64, λ::Float64)
    # If the bead coordinates are too far apart it's shorter to pass boundary
    for dim in eachindex(r′[begin, :])
        if norm(r′[begin, dim] - r′[end, dim]) > L
            # Shift r2 to nearest boundary of r1
            r′[end, dim] += sign.(r′[begin, dim] ) * (2 * L)
        end
    end

    # The path is constructed for 1+ particles
    m = size(r′, 1) - 2
    for j in 1:m
        # linear interpolation factor between the two boundary coordinates (Makes a nice curve)
        α = (m + 1 - j) / (m + 2 - j)

        # Construct a random displacement of a wiggle √(α*τ) and make sure to teleport it back within boundary
        r′[j+1, :] = α .* r′[j, :] + (1 - α) .* r′[end, :] + randn(size(r′, 2)) * √(2 * λ * α * τ)
    end
    for j in eachindex(@view r′[:, 1])
        r′[j, :] = teleport.(r′[j, :], L)
    end
end

function hardspherelevy!(r′::Matrix{Float64}, s::System, j₀::Int64, fpcycle::Function)::Bool
    # If the bead coordinates are too far apart it's shorter to pass boundary
    for dim in eachindex(r′[begin, :])
        if norm(r′[begin, dim] - r′[end, dim]) > s.L
            # Shift r2 to nearest boundary of r1
            r′[end, dim] += sign.(r′[begin, dim] ) * (2 * s.L)
        end
    end

    # The path is constructed for 1+ particles
    m = size(r′, 1) - 2
    for j in 1:m
        # linear interpolation factor between the two boundary coordinates (Makes a nice curve)
        α = (m + 1 - j) / (m + 2 - j)

        # Construct a random displacement of a wiggle √(α/2π)*λ_τ=√(2*λ*α*τ) and make sure to teleport it back within boundary. No new beads with a distance smaller than the 2d scattering length are created. We teleport the position only afterwards as the next position is didacted from the previous. We find the nn were we exclude the worldline that is altered
        pass = true; ctr = 0
        while pass
            pass = false; ctr +=1
            if ctr > s.ctr
                pass = true
                # @info "no new position outside hardsphere could be sample with levy"
                break
            end
            r′[j+1, :] = α .* r′[j, :] + (1 - α) .* r′[end, :] + randn(size(r′, 2)) * √(2 * s.λ * α * s.τ)
            # use the teleport position as these will be the final
            nn = find_nn(s, teleport.(r′[j+1, :], s.L), mod1(j₀ + j, s.M), exceptions = [fpcycle(mod1(j₀ + j, s.M))]) # j₀ + j because we start from j₀ + 1
            if nn != -1 && norm(distance.(teleport.(r′[j+1, :], s.L), s.world[nn].r[mod1(j₀ + j, s.M), :], s.L)) < s.a
                pass = true
            end
        end
        if pass
            return false
        end
    end
    # make sure every position is between periodic boundary
    for j in eachindex(@view r′[:, 1])
        r′[j, :] = teleport.(r′[j, :], s.L)
    end
    return true
end

function hardspherelevy!(r′::Matrix{Float64}, s::System, j₀::Int64)::Bool
    # If the bead coordinates are too far apart it's shorter to pass boundary
    for dim in eachindex(r′[begin, :])
        if norm(r′[begin, dim] - r′[end, dim]) > s.L
            # Shift r2 to nearest boundary of r1
            r′[end, dim] += sign.(r′[begin, dim] ) * (2 * s.L)
        end
    end

    # The path is constructed for 1+ particles
    m = size(r′, 1) - 2
    for j in 1:m
        # linear interpolation factor between the two boundary coordinates (Makes a nice curve)
        α = (m + 1 - j) / (m + 2 - j)

        # Construct a random displacement of a wiggle √(α/2π)*λ_τ=√(2*λ*α*τ) and make sure to teleport it back within boundary. No new beads with a distance smaller than the 2d scattering length are created. We teleport the position only afterwards as the next position is didacted from the previous. We find the nn were we exclude the worldline that is altered
        pass = true
        ctr = 0
        while pass
            ctr +=1; pass = false
            if ctr > s.ctr
                pass = true
                break
            end
            r′[j+1, :] = α .* r′[j, :] + (1 - α) .* r′[end, :] + randn(size(r′, 2)) * √(2 * s.λ * α * s.τ)
            nn = find_nn(s, teleport.(r′[j+1, :], s.L), mod1(j₀ + j, s.M)) # j₀ + j because we start from j₀ + 1
            if nn != -1 && norm(distance.(teleport.(r′[j+1, :], s.L), s.world[nn].r[mod1(j₀ + j, s.M), :], s.L)) < s.a
                pass = true
            end
        end
        if pass
            return false
        end
    end
    # make sure every position is between periodic boundary
    for j in eachindex(@view r′[:, 1])
        r′[j, :] = teleport.(r′[j, :], s.L)
    end
    return true
end

function sampleparticles(s::System, j₀::Int64, m::Int64)::Tuple{Int64,Int64}
    n1 = rand(1:s.N)
    if s.N==1
        error("we can only sample two particles if their are two particles")
    end
    # their must be the option taht n1=n2 so that it can get rejected based on their distance
    t = zeros(s.N)
    # compute the kinetic propagator for the beads (i, j₀) and (pcycle(j₀+m) of j, j₀+m)

    for i in 1:s.N
        Npol, pol = subcycle(s.world, i)
        nᵢnext = pcycle(j₀ + m, pol, Npol, s.M)
        # if s.world[nᵢnext] isa Worm
        #     t[i] = -Inf
        # else
        t[i] = s.lnK(
            s.world[n1].r[j₀, :],
            s.world[nᵢnext].r[mod1(j₀ + m, s.M), :],
            m * s.τ)
        # end
    end

    y = zeros(s.N)
    for i in 1:s.N
        Npol, pol = subcycle(s.world, n1)
        n1next =pcycle(j₀ + m, pol, Npol, s.M)
        # if s.world[n1next] isa Worm
        #     y[i] = -Inf
        # else
        y[i] = s.lnK(
            s.world[i].r[j₀, :],
            s.world[n1next].r[mod1(j₀ + m, s.M), :],
            m * s.τ)
        # end
    end

    norm = sum(exp.(t + y))

    particles = collect(1:s.N)
    weights = Weights(exp.(t + y) / norm)

    n2 = sample(particles, weights)
    return n1, n2
end

function prev(p::Vector{Worldline}, M::Int64, n₀::Int64, j₀::Int64)::Tuple{Int64,Int64}
    if iszero(n₀)
        return (0, 0)
    end
    n = j₀ == 1 ? cycle_findprev(p, n₀) : n₀
    j = mod1(j₀ - 1, M)
    if iszero(n)
        return (0, 0)
    elseif ismissing(p[n].r[j, 1])
        return (0, 0)
    else
        return (n, j)
    end
end
# function prev(::Vector{Worldline}, ::Int64, ::Missing, ::Missing)::Tuple{Missing,Missing}
    # (missing, missing)
# end

function next(p::Vector{Worldline}, M::Int64, n₀::Int64, j₀::Int64)::Tuple{Int64,Int64}
    if iszero(n₀)
        return (0, 0)
    end
    n = j₀ == M ? p[n₀].next : n₀
    j = mod1(j₀ + 1, M)
    if iszero(n)
        return (0, 0)
    elseif ismissing(p[n].r[j, 1])
        return (0, 0)
    else
        return (n, j)
    end
end

# function next(::Vector{Worldline}, ::Int64, ::Missing, ::Missing)::Tuple{Missing,Missing}
#     (missing, missing)
# end

function interaction_action!(w::Float64, s::System, fpcycle::Function, j::Int64)::Nothing
    nnlist = find_nns(s, fpcycle(j), mod1(j, s.M), exceptions=[fpcycle(j)])
    for nn in nnlist
        nn_nextbead = next(s.world, s.M, nn, mod1(j, s.M))
        j_nextbead = next(s.world, s.M, fpcycle(j), mod1(j, s.M))
        if iszero(nn_nextbead[1]) || iszero(j_nextbead[1])
            continue
        end
        r1 = distance.(
            s.world[nn].r[mod1(j, s.M), :],
            s.world[fpcycle(j)].r[mod1(j, s.M), :],
            s.L)
        r2 = distance.(
            s.world[nn_nextbead[1]].r[nn_nextbead[2], :],
            s.world[j_nextbead[1]].r[j_nextbead[2], :],
            s.L)
        w += s.lnU(r1, r2)
    end
end

function interaction_action!(w::Float64, s::System, r′::Matrix{Float64}, j::Int64, j₀::Int64, exceptions::Vector{Int64}=Int64[])::Nothing
    nnlist = find_nns(s, r′[j, :], mod1(j₀ + j - 1, s.M), exceptions=exceptions)
    for nn in nnlist
        nn_nextbead = next(s.world, s.M, nn, mod1(j₀ + j - 1, s.M))
        #TODO check if nn_nextbead[1] ==n1
        if iszero(nn_nextbead[1]) || nn_nextbead[1] ∈ exceptions
            continue
        end
        r1 = distance.(
            s.world[nn].r[mod1(j₀ + j - 1, s.M), :],
            r′[j, :],
            s.L)
        r2 = distance.(
            s.world[nn_nextbead[1]].r[nn_nextbead[2], :],
            r′[j+1, :],
            s.L)
            if 0.0 ∈ [r1...,r2...]
                println("r1=$r1 and r2=$r2")
            end
        w += s.lnU(r1, r2)
    end
end

function interaction_action!(w::Float64, s::System, r′::Array{Union{Missing, Float64}}, j::Int64, i::Int64, ii::Int64, mnext::Int64, exceptions::Vector{Int64}=Int64[])::Nothing
    nnlist = find_nns(s, disallowmissing(r′[j, :, i]), j, exceptions=exceptions)
    for nn in nnlist
        nn_nextbead = next(s.world, s.M, nn, j)
        if iszero(nn_nextbead[1])
            continue
        end
        r1 = distance.(
            s.world[nn].r[j, :],
            disallowmissing(r′[j, :, i]),
            s.L)
        r2 = distance.(
            s.world[nn_nextbead[1]].r[nn_nextbead[2], :],
            disallowmissing(r′[mnext, :, ii]),
            s.L)
        w += s.lnU(r1, r2)
    end
end

function move_polymer!(r′::Array{Union{Float64, Missing}}, s::System, maxd::Float64, pol::Vector{Int64})::Bool
    ctr = 0
    pass = true
    while pass
        pass = false
        ctr += 1
        if ctr > s.ctr
            return false
        end
        # Random displacement
        d = maxd * 2 .* (rand(s.dim) .- 0.5)
        # Shift the entire 'polymer' (all connected worldlines)
        for (i, n) in pairs(pol)
            for j in eachindex(skipmissing(s.world[n].r[:, 1]))
                r′[j, :, i] = teleport.(view(s.world[n].r, j, :) + d, s.L)

                # if the new distance is smaller than 2d scattering length, begin again
                nn = find_nn(s, disallowmissing(r′[j, :, i]), j, exceptions=[n])

                if nn != -1 && norm(distance.(r′[j, :, i], s.world[nn].r[j, :], s.L)) < s.a
                    pass = true
                    break # will only break the for loops (not the while)
                end
            end
        end
    end
    return true
end

# if I use only one update I get typeerror
# struct DummyUpdate <: UpdateC
#     counter::Counter
#     var::Step
#     counter_var::Counter

#     function DummyUpdate(s::System, step = 1.0; minstep = 1e-1, maxstep = s.L / 2, minacc = 0.4, maxacc = 0.6)
#         new(Counter(), Step(step, minstep, maxstep, minacc, maxacc), Counter())
#     end
# end

# function (u::DummyUpdate)(s::System)::Bool
#     false
# end

# mutable struct NumbOfLevels
#     l::Int64
#     minlevels::Int64
#     maxlevels::Int64
#     minacc::Float64
#     maxacc::Float64
# end


# function adjust!(m::NumbOfLevels, acceptance::Float64)
#     if acceptance < m.minacc
#         m.l -= 1
#     elseif acceptance > m.maxacc
#         m.l += 1
#     end

#     # Keep within bounds
#     m.l = max(m.minlevels, m.l)
#     m.l = min(m.maxlevels, m.l)
# end

# function expandbox(b1::Vector{Float64}, b2::Vector{Float64}, xold::Vector{Union{Missing,Float64}}, size::Float64)::Tuple{Vector{Float64},Vector{Float64}}
#     x10 = xold - b1
#     x02 = b2 - xold
#     x10 -= (2 * size) * round.(x10 / (2 * size))
#     x02 -= (2 * size) * round.(x02 / (2 * size))
#     return xold - x10, xold + x02
# end
# function expandbox(b1::Vector{Float64}, b2::Vector{Float64}, xold::Vector{Float64}, size::Float64)::Tuple{Vector{Float64},Vector{Float64}}
#     x10 = xold - b1
#     x02 = b2 - xold
#     x10 -= (2 * size) * round.(x10 / (2 * size))
#     x02 -= (2 * size) * round.(x02 / (2 * size))
#     return xold - x10, xold + x02
# end

# function sampleparticles(s::System, j₀::Int64, m::Int64, polW::Vector{Union{Missing,Int64}})::Tuple{Int64,Int64}
#     n1 = rand(1:s.N)
#     if n1 in skipmissing(polW)
#         return (-1, -1)
#     end
#    # their must be the option that n1=n2 so that it can get rejected based on their distance
#     t = Vector{Union{Float64,Missing}}(undef, s.N)
#     for j in 1:s.N
#         if j in skipmissing(polW)
#             t[j] = missing
#         else
#             Npol, pol = subcycle(s.world, j)
#             t[j] = s.lnK(
#                 view(s.world[n1].r, j₀, :),
#                 view(s.world[pcycle(j₀ + m, pol, Npol, s.M)].r, mod1(j₀ + m, s.M), :),
#                 m * s.τ)
#         end
#     end

#     y = zeros(s.N)
#     for j in 1:s.N
#         Npol, pol = subcycle(s.world, n1)
#         y[j] = s.lnK(
#             view(s.world[j].r, j₀, :),
#             view(s.world[pcycle(j₀ + m, pol, Npol, s.M)].r, mod1(j₀ + m, s.M), :),
#             m * s.τ)
#     end

#     norm = sum(skipmissing(exp.(t + y)))

#     particles = collect(eachindex(skipmissing(t)))
#     weights = Weights(disallowmissing(filter!(!ismissing, exp.(t + y) / norm)))

#     n2 = sample(particles, weights)
#     return n1, n2
# end
