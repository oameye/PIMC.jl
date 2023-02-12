using Pimc
using StatsBase

include(pwd()*"/examples/tools/logtools.jl")

potential = (r) -> 0.5*(r[1]^2+r[2]^2) # 2d harmonic oscillator
s = System(potential; dV=identity, λ=0.5, M=5, N=1, L=100.0, T=1.0, length_measurement_cycle=10);

n = 10_000; times = 10;
updates = [
    (1, SingleCenterOfMass(s, 1.0)),
    (1, ReshapeLinear(s, 2))
];
mea = ZMeasurement[Energy(s, 2*s.Ncycle*n*times)];

# Layout with general information for the initiation of the simulation
info(s)

# In the thermalization phase we take no measurements and we let the configuration converge
@info "thermalizing"
run!(s, n, updates)

while s.N_MC[s.N] < n*times
    run!(s, n, updates, Zmeasurements=mea)
end

# end # @everywhere
exact = 2.1639534137386534 # ground state energy
E_v = mean(collect(skipmissing(mea[1].energy_virial[s.N])))
E = mean(collect(skipmissing(mea[1].energy[s.N])))
