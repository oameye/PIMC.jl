module Pimc
using LinearAlgebra
using Distances: PeriodicEuclidean
using StatsBase
using Random
using Interpolations
using SpecialFunctions
using QuadGK
using NearestNeighbors
using Roots
using Optim
using DataStructures

# Main stuff
export System, run!, build_prop_int, worm_updates, acceptance, Coord, get_C

# Debug
export Worldline, Particle, Worm, levy!, distance, SectorCounter, update_nnbins!, check_nn, disallowmissing, debugpermlist, apply!, bin, lnV

# Updates
export permutations, subcycle, pcycle, Update, Updates
export Counter, Step, NumbOfLevels, NumbOfSlices, Threshold
export SingleBead, SingleCenterOfMass, PolymerCenterOfMass, ReshapeLinear, ReshapeSwapLinear, DummyUpdate, Update, UpdateWorm
export add_worm_updates, Insert, Remove, OpenGC, CloseGC, OpenHeadC, OpenTailC, CloseHeadC, CloseTailC, Advance, Recede, Swap, MoveHead, MoveTail, SwapC

#measurements
export Density, NumbOfParticle, Measurement, GMeasurement, ZMeasurement, Energy, ZG_ratio, SuperFluid

abstract type Worldline end
abstract type Measurement end
abstract type GMeasurement <: Measurement end
abstract type ZMeasurement <: Measurement end
abstract type Update end
abstract type UpdateWorm <: Update end
abstract type UpdateC <: Update end

allowmissing(x::AbstractArray{T}) where {T} = convert(AbstractArray{Union{T,Missing}}, x)
disallowmissing(x::AbstractArray{T}) where {T} = convert(AbstractArray{nonmissingtype(T)}, x)

include("stp.jl")
include("system.jl")
include("nn.jl")
include("updates/helper.jl")
include("updates/com.jl")
include("updates/reshape.jl")
include("updates/worms/worms_help.jl")
include("updates/worms/insert-remove.jl")
include("updates/worms/open-close.jl")
include("updates/worms/canonical/open-close_canonical.jl")
include("updates/worms/advance-recede.jl")
include("updates/worms/canonical/move_head-tail.jl")
include("updates/worms/swap_better.jl")
include("measurement.jl")

# Perform updates/measurements to the system.
# Each function comes with an integer that denotes the frequency of the Update.
Updates = Vector{Tuple{Int64,A}} where {A <: Update}
UpdatesWorm = Vector{Tuple{Int64,A}} where {A <: UpdateWorm}
GMeasurements = Vector{A} where {A <: GMeasurement}
ZMeasurements = Vector{A} where {A <: ZMeasurement}
Measurements =  Vector{A} where {A <: Measurement}

acceptance(q::Queue{Bool})::Float64 = sum(q)/length(q)

function queue!(c::Counter, acc::Bool)::Nothing
    c.tries += 1
    enqueue!(c.queue, acc ? 1 : 0)
    if length(c.queue) > c.range
        _ = dequeue!(c.queue)
    end
    nothing
end
function queue!(c::SectorCounter, acc::Bool)::Nothing
    enqueue!(c.queue, acc ? 1 : 0)
    if length(c.queue) > c.range
        _ = dequeue!(c.queue)
    end
    nothing
end

function measurement_Z_sector(s::System, measurements::ZMeasurements)::Nothing
    if isempty(measurements)
        return nothing
    end
    if s.worms == 0
        if any(isa.(s.world, Worm))
            error("s.worms = 0 but their isa worm")
        end
        if s.N ∉ keys(s.Nctr)
            s.Nctr[s.N]=0
            s.N_MC[s.N]=0
        end
        s.Nctr[s.N] += 1
        if s.Nctr[s.N] == s.Ncycle
            s.N_MC[s.N] += 1
            apply!(s, measurements)
            s.Nctr[s.N] = 0
        end
    end
    nothing
end

function measurement_G_sector(s::System, measurements::GMeasurements)::Nothing
    if isempty(measurements)
        return nothing
    end
    if 0 ∉ keys(s.Nctr)
        s.Nctr[0]=0
        s.N_MC[0]=0
    end
    s.Nctr[0] += 1
    if length(measurements)==1 && measurements[1] isa ZG_ratio
        apply!(s, measurements)
    elseif s.Nctr[0] % s.Ncycle == 0 && s.worms==0
        s.N_MC[0] += 1
        apply!(s, measurements)
    end
    nothing
end

function apply!(s::System, fs::Measurements)::Nothing
    for f in fs
        f(s)
    end
    nothing
end
function apply!(s::System, f::UpdateC)::Nothing
    acc = f(s)
    queue!(f.counter, acc)

    if f.counter_var.tries % f.counter_var.adj == 0
        adjust!(f.var, acceptance(f.counter_var.queue))
    end
    nothing
end
function apply!(s::System, f::UpdateWorm)::Nothing
    acc = f(s)
    queue!(f.counter, acc)

    if s.modifyW
        if f.counter_var.tries % f.counter_var.adj == 0
            adjust!(f.var, acceptance(f.counter_var.queue))
        end
    end
    nothing
end

# # Loop for updates and measurements.
# # In the thermalization phase we do not measurements
# function run!(s::System, n::Int64, updates::Updates)::Nothing
#     # number of tries to sample beads outside hardsphere
#     s.ctr = 10_000
#     (every, fs) = unzip(updates)
#     weights = Weights(1 ./ every)
#     for _ = 1:n
#         f = sample(fs, weights)
#         apply!(s, f)

#         # In thermelization phase we do not modify C (Z-G ratio) or count
#     end
#     nothing
# end

function run!(s::System, n::Int64, updates::Updates; Zmeasurements::ZMeasurements = ZMeasurement[], Gmeasurements::GMeasurements = GMeasurement[])::Nothing
    # number of tries to sample beads outside hardsphere; smaller than in thermalization phase
    therm = isempty(Zmeasurements) && isempty(Gmeasurements) ? true : false
    s.ctr = therm ? 10_000 : 1_000
    (every, fs) = unzip(updates)
    weights = Weights(1 ./ every)
    for _ = 1:n
        f = sample(fs, weights)
        apply!(s, f)

        if s.worms==0 && !s.gc && s.Ninit != s.N
            error("N=$(s.N)")
        end


        # equilibriated sampling between Z-G sector
        if s.worm_algorithm && !therm
            queue!(s.C, s.worms > 0)
            # if n % s.C.adj == 0 && s.C.modifyC && length(c.queue) > c.range
            #     adjust!(s.C)
            # end
        end

        measurement_G_sector(s, Gmeasurements)
        measurement_Z_sector(s, Zmeasurements)
    end
    nothing
end

function run!(s::System, n::Int64, zg_ratio::ZG_ratio)::Nothing
    # number of tries to sample beads outside hardsphere; smaller than in thermalization phase
    s.ctr = 1_000

    OC = Threshold(1, 1, s.M-1, 0.4, 0.6)
    # AR = Threshold(10, 1, s.M-1, 0.4, 0.6)
    worms = [
        (1,    Open(OC, adj = 1, range = 10_000)),
        (1,   Close(OC, adj = 1, range = 10_000)),
        # (1, MoveHead(AR, adj = 1, range = 10_000)),
        # (1, MoveTail(AR, adj = 1, range = 10_000))
    ]

    (every, fs) = unzip(worms)
    weights = Weights(1 ./ every)
    for _ = 1:n
        f = sample(fs, weights)
        _ = f(s)
        measurement_G_sector(s, [zg_ratio])
    end
    nothing
end

function unzip(a)
    (getfield.(a, x) for x in fieldnames(eltype(a)))
end

function add_worm_updates(s::System, updates::Updates;
    m::Int64 = 2,
    min::Int64 = 1, max::Int64 = s.M-1,
    minacc::Float64 = 0.4, maxacc::Float64 = 0.6,
    adj::Int64 = 1,
    range::Int64 = 10_000)::Tuple{Updates, Updates}
    if s.worm_algorithm
        IR = Threshold(m, min, max, minacc, maxacc)
        OC = Threshold(m, min, max, minacc, maxacc)
        AR = Threshold(m, min, max, minacc, maxacc)
        SW = NumbOfSlices(m, min, max, minacc, maxacc)
        if s.gc
            worms = [
                (1,  Insert(IR, adj = adj, range = range)),
                (1,  Remove(IR, adj = adj, range = range)),
                (1,    OpenGC(OC, adj = adj, range = range)),
                (1,   CloseGC(OC, adj = adj, range = range)),
                (1, Advance(AR, adj = adj, range = range)),
                (1,  Recede(AR, adj = adj, range = range)),
                (1,    SwapC(SW, adj = adj, range = range))
            ]
        else
            worms = [
                (1,    OpenGC(OC, adj = adj, range = range)),
                (1,   CloseGC(OC, adj = adj, range = range)),
                # (1,    OpenHeadC(OC, adj = adj, range = range)),
                # (1,   CloseHeadC(OC, adj = adj, range = range)),
                # (1,    MoveHead(AR, adj = adj, range = range))
                (1,    Insert(IR, adj = adj, range = range)),
                (1,    Remove(IR, adj = adj, range = range)),
                (1,   Advance(AR, adj = adj, range = range)),
                (1,    Recede(AR, adj = adj, range = range)),
                (1,     SwapC(SW, adj = adj, range = range))
            ]
        end
        updates = Vector{Tuple{Int64, Update}}(union(updates, worms))
    else
        error("Cannot add worm updates to a system which does not use worm algorithm")
    end
    return updates, worms
end

# function worm_updates(s::System;
#     m::Int64 = 5,
#     min::Int64 = 1, max::Int64 = s.M-1,
#     minacc::Float64 = 0.4, maxacc::Float64 = 0.6,
#     adj::Int64 = 1,
#     range::Int64 = 10_000)::Dict
#     if s.worm_algorithm
#         IR = Threshold(m, min, max, minacc, maxacc)
#         OC = Threshold(m, min, max, minacc, maxacc)
#         AR = Threshold(m, min, max, minacc, maxacc)
#         SW = NumbOfSlices(m, min, max, minacc, maxacc)
#         if s.gc
#             worms = Dict([
#                 (:insert,  Insert(IR, adj = adj, range = range)),
#                 (:remove,  Remove(IR, adj = adj, range = range)),
#                 (:open,    OpenGC(OC, adj = adj, range = range)),
#                 (:close,   CloseGC(OC, adj = adj, range = range)),
#                 (:advance, Advance(AR, adj = adj, range = range)),
#                 (:recede,  Recede(AR, adj = adj, range = range))
#                 # (1,    Swap(SW, adj = adj, range = range))
#             ])
#         else
#             worms = Dict([
#                 (:open,    OpenHeadC(OC, adj = adj, range = range)),
#                 (:close,   CloseHeadC(OC, adj = adj, range = range)),
#                 (:movehead, MoveHead(OC, adj = adj, range = range)),
#                 (:swap,     SwapC(SW, adj = adj, range = range))
#             ])
#         end
#     else
#         error("Cannot add worm updates to a system which does not use worm algorithm")
#     end
#     return worms
# end


# function run!(s::System, n::Int64, Zmeasurements::ZMeasurements, Gmeasurements::GMeasurements = GMeasurement[])::Nothing
#     s.ctr = 1_000
#     Zupdates = Dict([
#         (:com, PolymerCenterOfMass(s, 1.0)),
#         (:reshape, ReshapeLinear(s, 9))
#         # (:reshapeswap, ReshapeSwapLinear(s, 9))
#     ])
#     OC = Threshold(5, 1, s.M-1, 0.3, 0.7)
#     SW = NumbOfSlices(1, 1, s.M-1, 0.3, 0.7)
#     openclose = Dict([
#         (:open,    OpenHeadC(OC, adj = 1, range = 10_000)),
#         (:close,   CloseHeadC(OC, adj = 1, range = 10_000))

#     ])
#     Gupdates = Dict([
#         (:movehead, MoveHead(OC, adj = 1, range = 10_000)),
#         (:swap,     SwapC(SW, adj = 1, range = 10_000))
#     ])
#     Gupdates = worm_updates(s::System, m = 3)
#     for _ = 1:n
#         # apply!(s, rand(openclose)[2])
#         # if s.worm_algorithm
#         #     queue!(s.C, s.worms > 0)
#         # end
#         # if s.worms > 0
#         #     for _ in 1:5
#         #         apply!(s, rand(Gupdates)[2])
#         #     end
#         # end
#         if s.worms == 0
#             for _ in 1:15
#                 apply!(s, rand(Zupdates)[2])
#             end
#             measurement_G_sector(s, Gmeasurements)
#             measurement_Z_sector(s, Zmeasurements)
#         end
#     end
#     nothing
# end

# function debug(s::System)::Nothing
#     outofbound(s)
#     # worms_when_worms(s)
#     check_hard_sphere(s)
# end

# function check_nn(s::System, f::SubString{String})
#     for p in s.world
#         for j in eachindex(skipmissing(p.r[:, 1]))
#             if p.bins[j] != bin(disallowmissing(p.r[j, :]), s.nbins, s.L)
#                 error("Bin in wordline is not same as bin of position after function $(f).")
#             end
#         end
#     end
#     for j in eachindex(s.nn)
#         for b in eachindex(s.nn[j])
#             for n in s.nn[j][b]
#                 if !ismissing(s.world[n].r[j, 1])
#                     if bin(disallowmissing(s.world[n].r[j, :]), s.nbins, s.L) != b
#                         @info typeof(s.world[n])
#                         error("Bin of position $(bin(disallowmissing(s.world[n].r[j, :]), s.nbins, s.L)) is not in right bin of nn $(b) after function $(f)")
#                     end
#                     if s.world[n].bins[j] != b
#                         @info typeof(s.world[n])
#                         error("Bin in wordline is not in right bin of nn after function $(f)")
#                     end
#                 end
#             end
#         end
#     end
# end

# function check_nn(s::System)
#     for p in s.world
#         for j in eachindex(skipmissing(p.r[:, 1]))
#             if p.bins[j] != bin(disallowmissing(p.r[j, :]), s.nbins, s.L)
#                 error("Bin in wordline is not same as bin of position.")
#             end
#         end
#     end
#     for j in eachindex(s.nn)
#         for b in eachindex(s.nn[j])
#             for n in s.nn[j][b]
#                 if !ismissing(s.world[n].r[j, 1])
#                     if bin(disallowmissing(s.world[n].r[j, :]), s.nbins, s.L) != b
#                         @info typeof(s.world[n])
#                         error("Bin of position $(bin(disallowmissing(s.world[n].r[j, :]), s.nbins, s.L)) is not in right bin of nn $(b)")
#                     end
#                     if s.world[n].bins[j] != b
#                         @info typeof(s.world[n])
#                         error("Bin in wordline is not in right bin of nn")
#                     end
#                 end
#             end
#         end
#     end
# end

# function check_bead_link(s::System)::Nothing
#     for (n,p) in pairs(s.world)
#         for j in eachindex(skipmissing(p.V))
#             if ismissing(p.r[j, 1]) || ismissing(p.r[j, 2])
#                 error("the link vector of $n at time slice $j has a value but the position matrix has not")
#             end
#             if j == collect(eachindex(skipmissing(p.V)))[end]
#                 jnext = mod1(j+1, s.M)
#                 nnextbead = j == s.M ? p.next : n
#                 if ismissing(s.world[nnextbead].r[jnext, 1]) || ismissing(s.world[nnextbead].r[jnext, 2])
#                     error("the link vector of $n at time slice $j has a value but the position matrix has not")
#                 end
#             end
#         end
#         r_index = collect(eachindex(skipmissing(p.r[:,1])))
#         for j in r_index
#             if j == r_index[end]
#                 continue
#             else
#                 if ismissing(p.V[j]) || ismissing(p.r[j, 2])
#                     error("the link vector has a value but the position matrix has not")
#                 end
#             end

#         end
#     end
# end
# function check_bead_link(s::System, j₀, jₘ)::Nothing
#     for (n,p) in pairs(s.world)
#         for j in eachindex(skipmissing(p.V))
#             if ismissing(p.r[j, 1]) || ismissing(p.r[j, 2])
#                 error("j₀, jₘ = $(j₀),$(jₘ)
#                 \nthe link vector of $n at time slice $j has a value but the position matrix has not")
#             end
#             if j == collect(eachindex(skipmissing(p.V)))[end]
#                 jnext = mod1(j+1, s.M)
#                 nnextbead = j == 100 ? p.next : n
#                 if ismissing(s.world[nnextbead].r[jnext, 1]) || ismissing(s.world[nnextbead].r[jnext, 2])
#                     error("the link vector of $n at time slice $j has a value but the position matrix has not")
#                 end
#             end
#         end
#         r_index = collect(eachindex(skipmissing(p.r[:,1])))
#         for j in r_index
#             if j == r_index[end]
#                 continue
#             else
#                 if ismissing(p.V[j]) || ismissing(p.r[j, 2])
#                     error("the link vector has a value but the position matrix has not")
#                 end
#             end

#         end
#     end
# end


# function check_hard_sphere(s::System)::Nothing
#     for p in s.world
#         for q in s.world
#             for j in eachindex(view(p.r, :, 1))
#                 if p==q
#                     continue
#                 elseif ismissing(p.r[j,1]) || ismissing(q.r[j,1])
#                     continue
#                 elseif norm(distance.(p.r[j,:], q.r[j,:], s.L)) < s.a
#                     @error "Some particles are in eachothers hard sphere"
#                     error("Some particles are in eachothers hard sphere")
#                 end
#             end
#         end
#     end
# end

# function outofbound(s::System)::Nothing
#     for n in eachindex(s.world)
#         for j in s.M
#             for i in s.dim
#                 if !ismissing(s.world[n].r[j, i]) && ( s.world[n].r[j, i] <-s.L || s.L < s.world[n].r[j, i])
#                     @error "j=$(j) and n=$(n)"
#                 end
#             end
#         end
#     end
# end

# function worms_when_worms(s::System)::Nothing
#     if s.worms > 0
#         if s.worm_algorithm == false
#             @error "s.worms = $(s.worms); but s.worm_algorithm says: $(s.worm_algorithm)"
#         end
#     end
#     check=false
#     for p in Iterators.reverse(s.world)
#         if p isa Worm
#             check=true
#         end
#     end
#     if check && !s.worm_algorithm
#         @error "There is a worm: $check; but s.worm_algorithm says: $(s.worm_algorithm)"
#     end
#     if check && s.worms == 0
#         @error "There is a worm; but s.worms says: $(s.worms)"
#     end
#     if !check && s.worms > 0
#         @error "There is no worm; but s.worms says: $(s.worms)"
#     end
# end

function debugpermlist(s::System)
    perm = Vector{String}(undef, length(s.world))
    for (i,p) in pairs(s.world)
        if isa(p, Pimc.Worm)
            if !iszero(p.tail) && !iszero(p.head)
                perm[i] = "W"
            elseif iszero(p.tail)
                perm[i] = "H"
            elseif iszero(p.head)
                perm[i] = "T$(p.next)"
            else
                @error "Not a worm"
            end
        else
            perm[i] = "$(p.next)"
        end
    end
    return perm
end

# function wordline_zero(s::System, f::String)
#     for (i, p) in pairs(s.world)
#         if all(p.r .== 0.0)
#             error("particle $(i) is all zero after update $f")
#         end
#     end
# end

# function j₀100bug(s::System, f::String)
#     if s.worms == 0 || s.worm_algorithm == false
#         return nothing
#     end
#     ntail = rand(findall(i -> isa(i, Worm) && !iszero(i.tail), s.world))
#     NpolW, polW = subcycle(s.world, ntail)
#     nhead = polW[NpolW]

#     # take the lower and upper bound for the to-be modified time slices
#     j₀ = s.world[nhead].head
#     if j₀==100
#         @info "j₀=100 after update $f"
#     end
# end

end # module
