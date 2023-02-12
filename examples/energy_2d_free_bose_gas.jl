using Pimc, StatsBase

include(pwd()*"/examples/tools/logtools.jl")

N=1; M=10; T=1.0;

potential = _ -> 0.0
s = System(potential; dV=identity, λ=1.0, L=100.0, M=M, N=N, T=T,
           length_measurement_cycle=2, measure_scheme=:c)

n=20_000; times=10;
Zupdates = [
    (1, SingleCenterOfMass(s, 1.0)),
    (1, ReshapeLinear(s, 20))
    # (1, ReshapeSwapLinear(s, 20))
]
measurements =  ZMeasurement[
    Energy(s, s.Ncycle*n*times)
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
    run!(s, n, Zupdates; Zmeasurements = measurements)
end

@show mean(collect(skipmissing(measurements[1].energy[N])));
