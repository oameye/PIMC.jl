using Pimc
# using BenchmarkTools
# using JLD2, FileIO
# using Plots
# using StatsBase
# using OwnTime
# gr(show = true)
# Plots.GRBackend()

println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
# include(pwd()*"/examples/tools/plottools.jl")
# include("tools\\savetools.jl")
println("Ended precompiling")

function main()
    s = System(
        r -> 0
        # (r) -> 0.5*(r[1]^2+r[2]^2);
        dim = 2,
        M = 100,
        N = 2,
        L = 4.0,
        T = 1.0,
        worms = true,
        gc = false,
        measure_scheme = :gc,
        interactions = false)

    n = 10_000_000
    Zupdates = [
        (2, SingleCenterOfMass(s, 3.0)),
        (1, ReshapeLinear(s, 20))
    ]
    updates, Gupdates = add_worm_updates(s::System, Zupdates::Updates)

    measurents = [
        ZG_ratio()
    ]

    info(s)

    println("running")
    run!(s, n, Gupdates)

    # println("running")
    # run!(s, n, updates)

end
for _ in 1:10
    main()
end


# @code_warntype main()
# @profview main()
# owntime(stackframe_filter=filecontains("Pimc.jl"))


