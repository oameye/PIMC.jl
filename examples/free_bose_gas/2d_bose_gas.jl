using Pimc, StatsBase

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
        (1, SingleCenterOfMass(s, 1.0)),
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
