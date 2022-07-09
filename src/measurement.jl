VectorMissing{T} = Vector{Union{Missing,T}}

mutable struct Density <: ZMeasurement
    nbins::Int64 # amount of bins
    bin::Float64 # size of bin
    dens::Union{Array{Float64,1},Array{Float64,2}} #density matrix
    ndata::Int64 # Amount of density measurement taken

    Density(s::System; nbins=500) = new(nbins, (2 * s.L) / nbins, zeros(Float64, ntuple(i -> nbins, s.dim)), 0)
end

function (d::Density)(s::System)
    for p in s.world
        for m in 1:s.M
            ibin = floor.(Int64, (p.r[m, :] .+ s.L) ./ d.bin)
            if all(ibin .> 0) && all(ibin .< d.nbins + 1) # To avoid errors
                d.dens[Tuple(ibin)...] += 1
            end
        end
    end
    d.ndata += s.M
end

mutable struct NumbOfParticle <: GMeasurement
    particle::Vector{Union{Missing,Int64}}

    NumbOfParticle(n=20_000) = new(VectorMissing{Int64}(missing, n))
end

function (d::NumbOfParticle)(s::System)
    d.particle[findfirst(ismissing, d.particle)] = s.N
end

mutable struct Energy <: ZMeasurement
    energy::Dict{Int64,VectorMissing{Float64}}
    energy_virial::Dict{Int64,VectorMissing{Float64}}

    Energy(s::System, n=20_000) = new(
        Dict{Int64,VectorMissing{Float64}}([(0, VectorMissing{Float64}(missing, n))]),
        Dict{Int64,VectorMissing{Float64}}([(0, VectorMissing{Float64}(missing, n))])
    )
end

function (u::Energy)(s::System)
    p = s.world

    if s.worms > 0
        error("Energy can only be measured in the Z-sector.")
    end

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

mutable struct ZG_ratio <: GMeasurement
    Zsector::Int64
    total::Int64

    ZG_ratio() = new(0,0)
end

function (u::ZG_ratio)(s::System)
    u.total += 1
    if s.worms == 0
        u.Zsector += 1
    end
end

mutable struct SuperFluid <: ZMeasurement
    sf::Dict{Int64,VectorMissing{Float64}}
    # nbins::Int64 # amount of bins
    # bin::Float64 # size of bin
    # dens::Union{Array{Float64,1},Array{Float64,2}} #density matrix
    # ndata::Int64 # Amount of density measurement taken

    SuperFluid(s::System, n=20_000) = new(
        Dict{Int64,VectorMissing{Float64}}([(0, VectorMissing{Float64}(missing, n))])
        # Dict{Int64,VectorMissing{Float64}}([(0, VectorMissing{Float64}(missing, n))])
    )
end

function (u::SuperFluid)(s::System)
    p = s.world

    if s.worms > 0
        error("SuperFluidity can only be measured in the Z-sector.")
    end

    A, D = 0, 0
    for i in eachindex(p)
        for j in eachindex(p[i].r[:, 1])
            inext = j == s.M ? p[i].next : i
            jnext = mod1(j + 1, s.M)

            A = p[i].r[j, 1]*p[inext].r[jnext, 2] - p[i].r[j, 2]*p[inext].r[jnext, 1]
            D = p[i].r[j, 1]*p[inext].r[jnext, 1] + p[i].r[j, 2]*p[inext].r[jnext, 2]
        end
    end
    sf = (1/(s.λ*s.β))*(A/D)

    if s.N ∉ keys(u.sf)
        u.sf[s.N] = VectorMissing{Int64}(missing, length(u.sf[0]))
    end
    u.sf[s.N][findfirst(ismissing, u.sf[s.N])] = sf

end

# mutable struct VirialEnergy <: ZMeasurement
#     energy::Dict{Int64,Vector{Union{Missing,Float64}}}
#     # energy::Dict{Int64,Vector{Union{Missing,Float64}}}

#     VirialEnergy(s::System, n=20_000) = begin
#         if s.dV == zero
#             error("To measure the VirialEnergy you need the derivative of the external potential.")
#         end
#         new(
#         Dict{Int64,Vector{Union{Missing,Float64}}}([(0, Vector{Union{Missing,Int64}}(missing, n))])
#         # Dict{Int64,Vector{Union{Missing,Float64}}}([(0, Vector{Union{Missing,Int64}}(missing, n))])
#         )
#     end
# end

# function (u::VirialEnergy)(s::System)
#     p = s.world

#     if s.worms > 0
#         error("Energy can only be measured in the Z-sector.")
#     end

#     kin, pot = 0, 0
#     for i in eachindex(p)
#         for j in eachindex(p[i].r[:, 1])
#             inext = j == s.M ? p[i].next : i
#             jnext = mod1(j + 1, s.M)

#             kin += dot(p[i].r[j, :], s.dV(p[i].r[j, :]))
#             pot += s.V(p[i].r[j, :]) + s.V(p[inext].r[jnext, :])

#         end
#     end
#     E_s = 1 / (2 * s.M) * kin  + 1 / (2 * s.M) * pot

#     if s.N ∉ keys(u.energy)
#         u.energy[s.N] = Vector{Union{Missing,Int64}}(missing, length(u.energy[0]))
#     end
#     u.energy[s.N][findfirst(ismissing, u.energy[s.N])] = E_s


#      # pot, dpot, Qterm = 0, 0, 0, 0
#     # for i in eachindex(p)
#     #     inext =  p[i].next
#     #     iprev = cycle_findprev(p, i)

#     #     dr  =  p[i].r[s.M, :] - p[inext].r[1, :]

#     #     dr′ = [0,0]
#     #     for j in eachindex(p[i].r[:, 1])
#     #         if j != s.M
#     #             dr′ += sign.(p[i].r[j, :] -  p[i].r[j+1, :]) .* distance.(p[i].r[j, :], p[i].r[j+1, :], s.L)
#     #         else
#     #             dr′ += sign.(p[i].r[j, :] - p[inext].r[mod1(j+1, s.M), :]) .* distance.(p[i].r[j, :], p[inext].r[mod1(j+1, s.M), :], s.L)
#     #         end
#     #     end

#     #     Qterm +=  dot(dr, dr′)

#     #     for j in eachindex(p[i].r[:, 1])
#     #         if j != 1
#     #             dpot += dot(p[i].r[j, :]-p[i].r[1, :], s.dV(p[i].r[j, :]))
#     #         end

#     #         inext = j == s.M ? p[i].next : i
#     #         jnext = mod1(j + 1, s.M)

#     #         pot += 0.5*(s.V(p[i].r[j, :]) + s.V(p[inext].r[jnext, :]))
#     #     end
#     # end
#     # E = s.dim*s.N / (2 * s.β)
#     # E += 1 / (4 * s.λ * s.τ * s.M) * Qterm
#     # E += 1/ (2*s.M) * dpot
#     # E += pot

#     # if s.N ∉ keys(u.energy)
#     #     u.energy[s.N] = Vector{Union{Missing,Int64}}(missing, length(u.energy[0]))
#     # end
#     # u.energy[s.N][findfirst(ismissing, u.energy[s.N])] = E
# end

# function measure_energy(pimc::PIMC,par::MeasurementParams)
#     # Collect data from all time slices; All slices are equal :)
#     M = pimc.pimcpar.M
#     N = pimc.pimcpar.N
#     x = pimc.x
#     E  = 0
#     E2 = 0
#     for i in 1:M
#         for j in 1:N
#             xb = x[i,j]
#             # Virial theorem expression for energy
#             Ekin = 0.5*xb*dVdx(xb)
#             Epot = V(xb)
#             Evir = Ekin + Epot
#             E += Evir
#             E2 += Evir^2
#         end
#     end
#     par.data.E += E
#     par.data.E2 += E2
#     par.data.ndata += M
#     par
# end

# function get_energy(pimc::PIMC)
#     a = find_measurement(pimc,measure_energy)
#     p = a.par.data
#     E = p.E/p.ndata
#     E2 = p.E2/p.ndata
#     return E, E2
# end

# function save_energy(pimc::PIMC)
#     a = find_measurement(pimc,measure_energy)
#     p = a.par.data
#     E = p.E/p.ndata
#     E2 = p.E2/p.ndata
#     open(a.par.outfile,"a") do f
#         println(f, E," ",pimc.pimcpar.E_exact )
#     end
# end

# mutable struct Acceptance <: Measurement
#     particle::Vector{Union{Missing, Int64}}

#     Acceptance(s::System, updates::Vector{Update}) = new(Vector{Union{Missing, Int64}}(missing, n))
# end

# function (d::Acceptance)(s::System)
#     d.particle[findfirst(ismissing, d.particle)] = s.N
# end

#TODO radial distribution
#TODO Superfluid Fraction
#TODO Compressibilty