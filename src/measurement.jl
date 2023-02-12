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

    # if s.worms > 0
    #     error("Energy can only be measured in the Z-sector.")
    # end

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
