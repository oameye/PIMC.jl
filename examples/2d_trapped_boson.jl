# using Distributed
# using CSV
# using DataFrames
# addprocs(6)

# @everywhere begin
using Pimc
using StatsBase

# println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
include(pwd()*"/examples/tools/plottools.jl")
include(pwd()*"/examples/tools/savetools.jl")

# function main(M, T)
    potential = (r) -> 0.5*(r[1]^2+r[2]^2)
    s = System(potential;
                dV = identity,
                λ = 0.5,
                M = 5,
                N = 1, #!!!!!!!!!
                L = 100.0,
                T = 1.0,
                worms = false,
                length_measurement_cycle = 10
                );

    n = 10_000
    times = 10
    updates = [
        (1, SingleCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 2))
    ];
    mea = [
        Energy(s, 2*s.Ncycle*n*times)
        # Density(s)
    ];
    # Layout with general information for the initiation of the simulation
    info(s)
    # In the thermalization phase we take no measurements and we let the configuration converge
    println("thermalizing")
    run!(s, n, updates)

    while s.N_MC[s.N] < n*times
        run!(s, n, updates, Zmeasurements = mea)
    end
    # return (M,T, collect(skipmissing(measurements[1].energy[s.N])))
# end

# end # @everywhere
exact = 2.1639534137386534
E_v = mean(collect(skipmissing(mea[1].energy_virial[s.N])))
E = mean(collect(skipmissing(mea[1].energy[s.N])))
# plot_density(s, measurements[3])

# N=2
# z(b) = (exp(-0.5*b)/(1-exp(-b)))
# using ForwardDiff
# function energy(β)
    # -ForwardDiff.derivative(x -> log(z(x)^2), β)
# end;

# energy(1)
# 4.327906827477307


# Ms = [10, 50, 100, 200, 300, 400]
# Ts = [0.5, 1.0, 5.0, 10.0, 20.0]
# Ms = [10, 50]
# Ts = [0.5, 1.0]
# status = pmap(Iterators.product(Ms, Ts), on_error=identity) do (M, T)
#     main(M, T)
# end

# data = DataFrame()
# for out in status
#     data[!, "M$(out[1])T$(out[2])"] = out[3]
# end
# CSV.write("energy_2d_trapped_boson.csv", data)


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

# open("energy_2d_trapped_boson.txt", "w") do file
#         for out in status
#             println(file, out)
#         end
# end