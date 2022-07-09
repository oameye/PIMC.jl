using Distributed
using DataFrames
using CSV
str(x::Float64) = replace(string(x), "."=>"_")
addprocs(3)

@everywhere begin
using Pimc
using StatsBase
# println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
include(pwd()*"/examples/tools/plottools.jl")
include(pwd()*"/examples/tools/savetools.jl")

function main(M, N, T)
    potential = (r) -> 0.5*(r[1]^2+r[2]^2)
    s = System(potential;
                dV = identity,
                λ = 0.5,
                M = M,
                N = N,
                L = 100.0,
                T = T,
                measure_scheme = :c,
                length_measurement_cycle = 5
                )

    n = 20_000
    times = 50
    Zupdates = [
        (1, PolymerCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 20)),
        (1, ReshapeSwapLinear(s, 20))
    ]
    mea =  ZMeasurement[
        Density(s, nbins = 10000)
    ]
    # Layout with general information for the initiation of the simulation
    info(s)

    # In the thermalization phase we take no measurements and we let the configuration converge
    @info "thermalization"
    run!(s, times*n, Zupdates)

    @info "measurements"
    while s.N_MC[N] < times*n
        # info(updates)
        @info s.N_MC[N]
        run!(s, n, Zupdates; Zmeasurements = mea)
    end
    save_density(s, mea[1])
end
end # @everywhere
τ = 0.02
N = 5
pmap([0.05, 0.1, 0.5]) do T
    main(Int64(round(1/(τ*T))), N, T)
end

# τ = 0.05
# Ns = Int64[2, 3, 4, 5]
# Ts = Float64[5.0, 2.0, 1.0, 0.5, 0.1, 0.05, 0.01]
# E = pmap(Iterators.product(Ns, Ts), on_error=identity) do (N, T)
#     mea = main(Int64(round(1/(τ*T))), N, T);
#     Eᵥ = collect(skipmissing(mea[1].energy_virial[N]))
#     Eₜ = collect(skipmissing(mea[1].energy[N]))
#     l = length(Eᵥ)
#     matrix = hcat(ones(l) .* N, ones(l) .* T, Eᵥ, Eₜ)
#     return matrix
# end
# for i in eachindex(E)
#     N = E[i][1,1]
#     T = E[i][1,2]
#     df = DataFrame(E[i], ["N", "T", "energy_virial", "energy_thermo"])
#     CSV.write("trapped_bose_energy_N$(N)_T$(str(T)).csv", df)
# end

