# using Distributed
# addprocs(6)

# TODO Reshape worms: does nothing
# TODO Normal permutation updates: Get the right energies
# TODO Play with C (constant during whole sim): does not change anything
# TODO measure_scheme rigth?: seems ok

# @everywhere begin
using Pimc
using StatsBase

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

# println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
include(pwd()*"/examples/tools/plottools.jl")

function main(M, N, C)
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
                T = 1.0,
                SC = SC,
                worms = true,
                gc = false,
                modifyC = false,
                measure_scheme = :c,
                length_measurement_cycle = 10
                )

    n = 10_000
    times = 1
    Zupdates = [
        (1, PolymerCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 20))
        # (1, ReshapeSwapLinear(s, 20))
        #TODO make ReshapeSwapLinear Worm compatiable
    ]
    updates, Gupdates = add_worm_updates(s::System, Zupdates::Updates, m = 5)
    measurements =  ZMeasurement[
        VirialEnergy(s, s.Ncycle*n),
        Energy(s, s.Ncycle*n)
    ]
    # Layout with general information for the initiation of the simulation
    info(s)

    # In the thermalization phase we take no measurements and we let the configuration converge
    @info "thermalization"
    run!(s, times*n, Zupdates)
    run!(s, times*n, Gupdates)

    @info "measurements"
    # while s.N_MC[N] < n
    #     # info(updates)
    #     println(" ")
    #     @info get_C(s.C)
    #     @info acceptance(s.C.queue) #not N_g/N_z

    #     run!(s, times*n, updates, measurements, GMeasurement[])
    # end
    n_mc = s.N_MC[N]
    for i in 1:n
    run!(s, 1, updates, measurements, GMeasurement[])
        if n_mc != s.N_MC[N]
            N_E = findfirst(ismissing, measurements[2].energy[N]) -1
            E = measurements[2].energy[N][N_E]
            if E < 0
                println("negative: $i\t energy: $E")
                n_mc = s.N_MC[N]
                display(worldlines(s))
                wait_for_key("press a key")
            end
            if E > 0
                println("positive: $i\t energy: $E")
                n_mc = s.N_MC[N]
                display(worldlines(s))
                wait_for_key("press a key")
            end
        end
    end

    return measurements
end

# [det_C(0.2, 10, 1.0) for _ in 1:5]
mea = main(10, 2, 3.0);
E_v = mean(skipmissing(mea[1].energy[2]))
E = mean(skipmissing(mea[2].energy[2]))
# E_ex = 4.028429828781551
# end # @everywhere


# Ms = [10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
# Ms = [10, 50]
# status = pmap(Ms, on_error=identity) do (M)
    # main(M)
# end
# status[1][1]

# open("energy_2d_trapped_boson.txt", "w") do file
#     header = "N_MC, "
#     for out in status
#         header *= "M$(out[1]), "
#     end
#     println(file, header)
#     for i in eachindex(status[1][2])
#         line ="$i, "
#         for out in status
#             line *= "$(out[2][i]), "
#         end
#         println(file, line)
#     end
# end
# include("plottools.jl")
# include("savetools.jl")
# function det_C(C,M,T)
#     potential = (r) -> 0.5*(r[1]^2+r[2]^2)
#     s = System(
#     # potential;
#     r -> 0,
#     dim = 2,
#     M = M,
#     N = 1,
#     L = 100.0,
#     T = T,
#     C = C,
#     worms = true,
#     gc = false,
#     measure_scheme = :gc,
#     length_measurement_cycle=5,
#     modifyC=false)

#     n = 1_000_000

#     ratio = ZG_ratio()

#     run!(s, n, ratio)
#     # println("C = $(Pimc.C(s.C))")
#     # println("total = $(ratio.total)")
#     # println("Zsector = $(ratio.Zsector)")
#     return (ratio.total-ratio.Zsector)/ratio.Zsector
# end