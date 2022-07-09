using Pimc
using StatsBase
using Plots

# println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
include(pwd()*"/examples/tools/plottools.jl")

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

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
                L = 4.0,
                T = 1.0,
                SC = SC,
                worms = true,
                gc = false,
                modifyC = false,
                measure_scheme = :c,
                length_measurement_cycle = 1
                )

    n = 1_000
    times = 3
    Zupdates = [
        (1, PolymerCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 20)),
        (1, ReshapeSwapLinear(s, 20))
    ]
    updates, Gupdates = add_worm_updates(s::System, Zupdates::Updates, m = 5)
    # mea =  ZMeasurement[
    #     VirialEnergy(s, s.Ncycle*n),
    #     Energy(s, s.Ncycle*n)
    # ]

    # Layout with general information for the initiation of the simulation
    info(s)

    # In the thermalization phase we take no measurements and we let the configuration converge
    @info "thermalization"
    run!(s, times*n, Zupdates)
    # run!(s, times*n, Gupdates)

    # function nextframe(s, updates)
    #     run!(s, 1, updates)
    #     return worldlines(s)
    # end
    # run!(s, 100, Gupdates, mea)
    for _ in 1:n
        run!(s, 1, Zupdates)
        # E =  mea[2].energy[N][findfirst(ismissing, mea[2].energy[N])-1]
        # if E < 0.0
        #     println("E=$E")
        #     display(worldlines(s))
        #     wait_for_key("wait_for_key")
        # end
        display(topview(s))
        wait_for_key("wait_for_key")
    end
    # anim = @animate for _ in 1:n
    #     nextframe(s, updates)
    # end
    # mp4(anim, "test.mp4", fps = 50)
    # display(worldlines(s))
end

# end # @everywhere
mea = main(100, 2, 0.01);




