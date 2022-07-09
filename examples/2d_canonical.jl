using Distributed
addprocs(6)


@everywhere begin
using Pimc
# using StatsBase
println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl");
include(pwd()*"/examples/tools/potentialtools.jl");
include(pwd()*"/examples/tools/plottools.jl");
include(pwd()*"/examples/tools/savetools.jl");
println("Ended precompiling")

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

function main(L, N, g, M, T, n, times, potential)
    propint = build_prop_int(round(√(L^2+L^2), RoundUp), g, 1/(T*M))
    s = System(
        potential;
        λ = 1/(π^2),
        M = M,
        N = N,
        L = L,
        T = T,
        propint = propint,
        interactions = true,
        worms = false,
        gc = false,
        length_measurement_cycle = 3
        )

    updates = [
        (1, SingleCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 5)),
        (1, ReshapeSwapLinear(s, 20))
    ]
    mea = [
        Density(s)
    ]

    info(s) # Layout with general information for the initiation of the simulation

    # In the thermalization phase we take no measurements and we let the configuration converge
    println("Thermalizing")
    run!(s, n*times, updates)

    println("Running")
    while s.N_MC[s.N] < n*times
        run!(s, n, updates, Zmeasurements=mea)
        @info s.N_MC[s.N]
    end

    return s, mea[]
end

# cyclic_4 = [2π * k / 4 for k in 0:3];
# x = collect(range(-8,8, length=500))
# y = collect(range(-8,8, length=500))
# X, Y = meshgrid(x, y)
# V = Matrix{Float64}(undef, 500,500)
# for i in eachindex(V)
#     V[i] = potential([X[i], Y[i]])
# end
# heatmap(x, y, V)

end

gs = [2.0, 20.0]
depths = [4.0 ,6.0, 10.0]
pmap(Iterators.product(gs, depths)) do (g, V₀)
    scale = 1.0;
    name = :l25
    potential = generate_V(scale, V₀, name; attractive = true)
    n = 100_000; times = 10
    L = 8.0
    M = 200
    N = 20
    T = 0.2
    s, d = main(L, N, g, M, T, n, times, potential);
    save_density(s, d, g, name, V₀)
end



# plot_density(s, mea[1])
# @time main()
# @benchmark main()
# @code_warntype main()
# @profview main()
# owntime(stackframe_filter=filecontains("Pimc.jl"))