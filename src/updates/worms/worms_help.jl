mutable struct Threshold
    m::Int64
    min::Int64
    max::Int64
    minacc::Float64
    maxacc::Float64
end

function adjust!(m::Threshold, acc::Float64)
    if acc < m.minacc
        m.m -= 1 #Notice, inverse of NumbOfSlices!
    elseif acc > m.maxacc
        m.m += 1
    end

    # Keep within bounds
    m.m = max(m.min, m.m)
    m.m = min(m.max, m.m)
end

function lengthpolymer(p::Vector{Worldline}, pol::Vector{Int64}, sym::Symbol=:bead)::Int64
    m = 0
    if sym == :bead
        for n in pol
            m += count(.!isa.(view(p[n].r, :, 1), Missing))
        end
    else sym == :link
        for n in pol
            m += count(.!isa.(p[n].V, Missing))
        end
    end
    return m
end

function fixindexing!(p::Vector{Worldline}, n::Vector{Int64})::Nothing
    for l in p
        if !iszero(l.next)
            l.next -= count(<(l.next), n)
        end
    end
end

function fixindexing!(p::Vector{Worldline}, nn::NearestNeighbors, n::Vector{Int64})::Nothing
    for (i,l) in pairs(p)
        if !iszero(l.next)
            l.next -= count(<(l.next), n)
        end
        if i ∉ n
            cc = count(<(i), n)
            if cc > 0
                fixindexing_nn!(nn, i, i-cc, l.bins)
            end
        end
    end
end

# function sampleparticle_swap(s::System, m::Int64, n₀::Int64, j₀::Int64, p_nnbins::Vector{Int64})::Tuple{Int64,Float64}
#     t = createKlist(s, m, n₀, j₀, p_nnbins)
#     norm = sum(t)
#     t /= norm

#     n = sample(p_nnbins, Weights(t))
#     return n, norm
# end

# function createKlist(s::System, m::Int64, n₀::Int64, j₀::Int64, p_nnbins::Vector{Int64})::Vector{Float64}
#     t = Vector{Float64}(undef, length(p_nnbins))
#     # do not need to check for ismissing as it is already done for the bins
#     for (i, p) in pairs(p_nnbins)
#         t[i] = exp(
#             s.lnK(
#                 disallowmissing(s.world[n₀].r[j₀, :]),
#                 disallowmissing(s.world[p].r[mod1(j₀ + m, s.M), :]),
#                 m * s.τ)
#         )
#     end
#     return t
# end

function sampleparticle_swap(s::System, m::Int64, n₀::Int64, j₀::Int64)::Tuple{Int64,Float64}
    t = createKlist(s, m, n₀, j₀)

    norm = 0.0
    for i in t
        norm += i
    end

    for i in eachindex(t)
        t[i] /= norm
    end
    particles = collect(eachindex(t))
    weights = Weights(t)

    n = sample(particles, weights)
    return n, norm
end

function createKlist(s::System, m::Int64, n₀::Int64, j₀::Int64)::Vector{Float64}
    t = zeros(Float64, s.N)
    for i in eachindex(t)
        if !ismissing(s.world[i].r[mod1(j₀ + m, s.M), 1])
            t[i] = exp(
                s.lnK(
                    disallowmissing(s.world[n₀].r[j₀, :]),
                    s.world[i].r[mod1(j₀ + m, s.M), :],
                    m * s.τ
                )
            )
        end
    end
    return t
end

function gaussian_sample_hardsphere(s::System, r₀::Vector{Float64}, j::Int64, m::Int64; exceptions::Vector{Int64}  = Int64[])::Vector{Float64}
    pass = true; ctr = 0
    while pass
        pass = false; ctr += 1
        if ctr > s.ctr
            pass = true
            # @info "no new position outside hardsphere could be sampled with gaussian distribution"
            break
        end

        r′ =  r₀ + randn(s.dim) .* √(2 * s.λ * m * s.τ)
        # use the teleport position as these will be the final
        nn = find_nn(s, teleport.(r′, s.L), mod1(j, s.M), exceptions=exceptions)
        if nn != -1 && norm(distance.(teleport.(r′, s.L), s.world[nn].r[mod1(j, s.M), :], s.L)) < s.a
            pass = true
        end
    end

    if pass
        r′ = zeros(Float64, 2)
    end
    return r′
end

function uniform_sample_hardsphere(s::System, j::Int64, exceptions::Vector{Int64} = Int64[])::Vector{Float64}
    pass = true; ctr = 0
    while pass
        pass = false; ctr += 1
        if ctr > s.ctr
            pass = true
            break
        end
        r = (2 * s.L) .* (rand(Float64, s.dim) .- 0.5)
        # use the teleport position as these will be the final
        nn = find_nn(s, teleport.(r, s.L), mod1(j, s.M), exceptions=exceptions)
        if nn != -1 && norm(distance.(teleport.(r, s.L), s.world[nn].r[mod1(j, s.M), :], s.L)) < s.a
            pass = true
        end
    end

    if pass
        r = zeros(Float64, 2)
    end

    return r
end