# using Distributed
# addprocs(5)

# TODO Reshape worms: does nothing
# TODO Normal permutation updates: Get the right energies
# TODO Play with C (constant during whole sim): does not change anything
# TODO measure_scheme rigth?: seems ok

# FIXME ReshapeSwapLinear with worms?
# FIXME one boson worms without swap, right energy? otherwise something worng with worms without swap

# @everywhere begin
using Pimc
using StatsBase
# using Plots

# println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")

function main(M, T, N, C)
    SC = SectorCounter(C,
        min = 1e-6, max = 1e1,
        range = 20_000,
        minacc = 0.30, maxacc = 0.70,
        modifyC = false
        )

    potential = (r) -> 0.5*(r[1]^2+r[2]^2)
    s = System(potential;
                dV = identity,
                λ = 0.5,
                M = M,
                N = N,
                L = 100.0,
                T = T,
                SC = SC,
                worms = true,
                gc = false,
                modifyC = false,
                measure_scheme = :c,
                length_measurement_cycle = 1
                );

    n = 10_000
    times = 10
    Zupdates = [
        (1, PolymerCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 5)),
        (1, ReshapeSwapLinear(s, 5))
    ];
    updates, Gupdates = add_worm_updates(s::System, Zupdates::Updates, m = 5);
    measurements =  ZMeasurement[
        VirialEnergy(s, 2*s.Ncycle*n*times),
        Energy(s, 2*s.Ncycle*n*times)
    ];
    # Layout with general information for the initiation of the simulation
    info(s)

    # In the thermalization phase we take no measurements and we let the configuration converge
    @info "thermalization"
    run!(s, times*n, Zupdates)

    @info "measurements"
    while s.N_MC[N] < n*times
        # info(updates)
        println(" ")
        @info s.N_MC[N]
        @info get_C(s.C)
        @info acceptance(s.C.queue) #not N_g/N_z

        run!(s, n, measurements)
    end
    return measurements
end

# end # @everywhere
M = 10
N = 2
T = 1.0
C = 0.5
mea = main(M, T, N, C);
Eᵥ = mean(skipmissing(mea[1].energy[N]))
E = mean(skipmissing(mea[2].energy[N]))
println("Eᵥ = $(Eᵥ)")
println("E = $(E)")

# status = pmap(zeros(5)) do _
#     mea = main(10, 1.0, 2e-4);
#     Eᵥ = mean(skipmissing(mea[1].energy[1]))
#     E = mean(skipmissing(mea[2].energy[1]))
#     return [Eᵥ, E]
# end
# E = Matrix(reduce(hcat,status)')
# println(mean(E[:, 1]))
# println(mean(E[:, 2]))
# exact = 2.1639534137386534
# 2.1605948385183447 worms open-close (long run) ✔
# 2.160428482276161 only canonical ✔
# 2.168085913795338 worms open-close-MoveHead-MoveTail ✖
