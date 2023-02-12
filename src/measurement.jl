function measurement_Z_sector(s::System, measurements::ZMeasurements)::Nothing
    if isempty(measurements)
        return nothing
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
    nothing
end

"""

The `Density` struct is a mutable struct that contains information about a density measurement.
It is a subtype of the `ZMeasurement` struct. The struct contains the following fields:
    - `nbins`: an `Int64` representing the amount of bins in the density measurement.
    - `bin`: a `Float64` representing the size of each bin in the density measurement.
    - `dens`: a `Union` of either an `Array{Float64, 1}` or an `Array{Float64, 2}` representing the density matrix. The dimensions of the array will depend on the number of dimensions of the system being measured (`s.dim`).
    - `ndata`: an `Int64` representing the amount of density measurements taken.

    The `Density` struct also contains a constructor that takes in a `System` object and an optional keyword argument `nbins`. The constructor creates a new `Density` object with `nbins` as the number of bins, `(2 * s.L) / nbins` as the size of each bin, a `dens` array of zeros with dimensions determined by `ntuple(i -> nbins, s.dim)`, and `ndata` set to zero.
"""

mutable struct Density <: ZMeasurement
    nbins::Int64 # amount of bins
    bin::Float64 # size of bin
    dens::Union{Array{Float64,1},Array{Float64,2}} # density matrix
    ndata::Int64 # Amount of density measurement taken

    Density(s::System; nbins=500) = new(nbins, (2 * s.L) / nbins, zeros(Float64, ntuple(i -> nbins, s.dim)), 0)
end

"""
The function takes in a `Density` object (`d`) and a `System` object (`s`) as input. It loops over all particles in the `s.world` and calculates the bin index (`ibin`) for each particle. If the bin index for all dimensions is greater than 0 and less than `d.nbins + 1`, the density at that bin is incremented by 1. The `ndata` field is then incremented by `s.M` to keep track of the total number of density measurements taken.

"""

function (d::Density)(s::System)
    for p in s.world
        for m in 1:s.M
            ibin = floor.(Int64, (p.r[m, :] .+ s.L) ./ d.bin)
            if all(ibin .> 0) && all(ibin .< d.nbins + 1)
                d.dens[Tuple(ibin)...] += 1
            end
        end
    end
    d.ndata += s.M
end


"""
`Energy` is a mutable structure which is a subtype of `ZMeasurement`. It contains two dictionaries `energy` and `energy_virial`.

### Fields

- `energy::Dict{Int64,VectorMissing{Float64}}`: A dictionary where the key is the number of polymer chains and the value is a `VectorMissing{Float64}` of size n.
- `energy_virial::Dict{Int64,VectorMissing{Float64}}`: A dictionary where the key is the number of polymer chains and the value is a `VectorMissing{Float64}` of size n.

### Constructor

```julia
Energy(s::System, n=20_000) = new(
    Dict{Int64,VectorMissing{Float64}}([(0, VectorMissing{Float64}(missing, n))]),
    Dict{Int64,VectorMissing{Float64}}([(0, VectorMissing{Float64}(missing, n))])
)
```

The constructor takes in a `System` object and an optional argument `n` which sets the size of the `Vector{Union{Missing,T}}` objects in the dictionaries to `n`. The constructor initializes the two dictionaries with the key `0` and a `Vector{Union{Missing,T}}` object of size n with all values as `missing`.

"""
mutable struct Energy <: ZMeasurement
    energy::Dict{Int64,VectorMissing{Float64}}
    energy_virial::Dict{Int64,VectorMissing{Float64}}

    Energy(s::System, n=20_000) = new(
        Dict{Int64,VectorMissing{Float64}}([(0, VectorMissing{Float64}(missing, n))]),
        Dict{Int64,VectorMissing{Float64}}([(0, VectorMissing{Float64}(missing, n))])
    )
end

"""
    This function calculates the energy and virial energy of a given `System` object and updates the `Energy` type. The function first defines four local variables `link`, `pot`, `vkin`, and `vpot` to store intermediate values. Then, it loops over the particles in the system and calculates the energy and virial energy based on the interactions between particles. The final values of energy and virial energy are stored in the `u.energy` and `u.energy_virial` fields of the `Energy` object, respectively. If the number of particles in the system `s.N` is not present in the energy dictionary, a new key is added and the energy value is added to the corresponding vector.
"""

function (u::Energy)(s::System)
    p = s.world

    link, pot = 0, 0
    vkin, vpot = 0, 0
    for i in eachindex(p)
        for j in eachindex(p[i].r[:, 1])
            inext = j == s.M ? p[i].next : i
            jnext = mod1(j + 1, s.M)

            dr = distance.(p[i].r[j, :], p[inext].r[jnext, :], s.L)

            link += dot(dr, dr)
            pot += s.V(p[i].r[j, :]) + s.V(p[inext].r[jnext, :])

            vkin += dot(p[i].r[j, :], s.dV(p[i].r[j, :]))
            vpot += s.V(p[i].r[j, :]) + s.V(p[inext].r[jnext, :])
        end
    end
    E = s.dim*s.N / (2 * s.τ) - 1 / (4 * s.λ * s.τ^2 * s.M) * link + 1 / (2 * s.M) * pot
    E_v = 1 / (2 * s.M) * vkin  + 1 / (2 * s.M) * vpot

    if s.N ∉ keys(u.energy)
        u.energy[s.N] = VectorMissing{Int64}(missing, length(u.energy[0]))
    end
    if s.N ∉ keys(u.energy_virial)
        u.energy_virial[s.N] = VectorMissing{Int64}(missing, length(u.energy_virial[0]))
    end
    u.energy[s.N][findfirst(ismissing, u.energy[s.N])] = E
    u.energy_virial[s.N][findfirst(ismissing, u.energy_virial[s.N])] = E_v
end


#TODO radial distribution
#TODO Superfluid Fraction
#TODO Compressibilty
