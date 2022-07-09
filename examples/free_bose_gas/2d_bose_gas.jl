# FIXME if I increase L I get the right result. If I decrease L I get the wrong energy levels, why?

using Pimc
using StatsBase
# println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")

function main(M, N, T)
    potential = _ -> 0.0
    s = System(potential;
                dV = identity,
                λ = 1.0,
                M = M,
                N = N,
                L = 100.0,
                T = T,
                measure_scheme = :c,
                length_measurement_cycle = 2
                )

    n = 20_000
    times = 10
    Zupdates = [
        (1, PolymerCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 20))
        # (1, ReshapeSwapLinear(s, 20))
    ]
    mea =  ZMeasurement[
        Energy(s, s.Ncycle*n*times)
        # SuperFluid(s, 2*times*n*s.Ncycle)
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
    return mea
end
# end # @everywhere
M = 10; N = 1; T = 1.0
mea = main(M, N, T);
@show mean(collect(skipmissing(mea[1].energy[N])));

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

