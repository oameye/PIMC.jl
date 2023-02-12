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

#measurements
export Density, NumbOfParticle, Measurement, GMeasurement, ZMeasurement, Energy, ZG_ratio, SuperFluid

abstract type Worldline end
abstract type Measurement end
abstract type ZMeasurement <: Measurement end
abstract type Update end
abstract type UpdateC <: Update end

allowmissing(x::AbstractArray{T}) where {T} = convert(AbstractArray{Union{T,Missing}}, x)
disallowmissing(x::AbstractArray{T}) where {T} = convert(AbstractArray{nonmissingtype(T)}, x)

include("propagator.jl")
include("system.jl")
include("nearest_neighbours.jl")
include("updates/helper.jl")
include("updates/com.jl")
include("updates/reshape.jl")
include("measurement.jl")

# Perform updates/measurements to the system.
# Each function comes with an integer that denotes the frequency of the Update.
Updates = Vector{Tuple{Int64,A}} where {A <: Update}
ZMeasurements = Vector{ZMeasurement}
Measurements = Vector{Measurement}

acceptance(q::Queue{Bool})::Float64 = sum(q)/length(q)

function queue!(c::Counter, acc::Bool)::Nothing
    c.tries += 1
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
    # if s.worms == 0
    # if any(isa.(s.world, Worm))
    #     error("s.worms = 0 but their isa worm")
    # end
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
    # end
    nothing
end

function apply!(s::System, fs::Vector{Measurement})::Nothing
    for f in fs
        f(s)
    end
    nothing
end
apply!(s::System, fs::ZMeasurements) = apply!(s,Measurements(fs))
function apply!(s::System, f::UpdateC)::Nothing
    acc = f(s)
    queue!(f.counter, acc)

    if f.counter_var.tries % f.counter_var.adj == 0
        adjust!(f.var, acceptance(f.counter_var.queue))
    end
    nothing
end

function run!(s::System, n::Int64, updates::Updates; Zmeasurements::ZMeasurements = ZMeasurement[])::Nothing
    # number of tries to sample beads outside hardsphere; smaller than in thermalization phase
    therm = isempty(Zmeasurements) ? true : false
    s.ctr = therm ? 10_000 : 1_000
    (every, fs) = unzip(updates)
    weights = Weights(1 ./ every)
    for _ = 1:n
        f = sample(fs, weights)
        apply!(s, f)

        measurement_Z_sector(s, Zmeasurements)
    end
    nothing
end

function unzip(a)
    (getfield.(a, x) for x in fieldnames(eltype(a)))
end

end # module
