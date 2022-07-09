# using Distributed
# addprocs(10)


# @everywhere begin
using Pimc
using StatsBase
# using Plots

# println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")

function main(M, T, C)
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
                N = 1,
                L = 100.0,
                T = T,
                SC = SC,
                worms = true,
                gc = false,
                modifyC = false,
                measure_scheme = :c,
                length_measurement_cycle = 3
                )

    n = 20_000
    times = 5
    Zupdates = [
        (1, SingleCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 5))
        # (1, ReshapeSwapLinear(s, 20))
    ]
    updates, Gupdates = add_worm_updates(s::System, Zupdates::Updates, m = 5)
    mea =  ZMeasurement[
        Energy(s, s.Ncycle*n*times)
    ]
    # Layout with general information for the initiation of the simulation
    info(s)

    # In the thermalization phase we take no measurements and we let the configuration converge
    @info "thermalization"
    run!(s, times*n, Zupdates)
    run!(s, times*n, updates)
    # run!(s, times*n, Gupdates)

    @info "measurements"
    while s.N_MC[1] < n*times
        # info(updates)
        println(" ")
        @info s.N_MC[1]
        @info get_C(s.C)
        @info acceptance(s.C.queue) #not N_g/N_z

        run!(s, n, updates; Zmeasurements = mea)
    end
    return mea
end

# end # @everywhere
mea = main(10, 1.0, 4e-2);
Eᵥ = mean(skipmissing(mea[1].energy_virial[1]))
E = mean(skipmissing(mea[1].energy[1]))
# status = pmap(zeros(10)) do _
#     mea = main(20, 1.0, 5e-5);
#     Eᵥ = mean(skipmissing(mea[1].energy[1]))
#     E = mean(skipmissing(mea[2].energy[1]))
#     return [Eᵥ, E]
# end
# E = Matrix(reduce(hcat,status)')
# mean(E[:, 1])
# mean(E[:, 2])
# exact = 2.1639534137386534
# 2.1605948385183447 worms open-close (long run) ✔
# 2.160428482276161 only canonical ✔
# 2.168085913795338 worms open-close-MoveHead-MoveTail ✖

# if !isa(s.world[nhead], Worm)
#     println(pol)
#     println(s.N)
#     println("$(debugpermlist(s))")
#     println("debug")
# end

