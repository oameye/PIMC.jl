# using Distributed
# addprocs(22)

# @everywhere begin
using Pimc
using StatsBase

using DataFrames, CSV
str(x::Float64) = replace(string(x), "."=>"_")

# println("Compiling plot and save function")includ
include(pwd()*"/examples/tools/logtools.jl");
# include(pwd()**"/examples/tools/energytools.jl")

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
        length_measurement_cycle = 3
        )

    n = 200_000
    times = 5
    Zupdates = [
        (1, PolymerCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 5)),
        (1, ReshapeSwapLinear(s, 5))
    ]
    updates, Gupdates = add_worm_updates(s::System, Zupdates::Updates, m = 3)
    mea =  ZMeasurement[
        Energy(s, s.Ncycle*times*n)
    ]
    # Layout with general information for the initiation of the simulation
    info(s)

    # In the thermalization phase we take no measurements and we let the configuration converge
    @info "thermalization"
    run!(s, times*n, Zupdates)
    # run!(s, times*n, Gupdates)

    @info "measurements"
    while s.N_MC[N] < times*n
        # info(updates)
        println(" ")
        @info s.N_MC[N]
        @info get_C(s.C)
        @info acceptance(s.C.queue) #not N_g/N_z

        run!(s, n, updates; Zmeasurements = mea)
    end
    return mea
end
# end # @everywhere
M = 20
T = 1.0
N = 10
C = 3e-5
mea = main(M, T, N, C);
@show Eᵥ = mean(skipmissing(mea[1].energy_virial[N]));
@show δEᵥ = std(skipmissing(mea[1].energy_virial[N]));
@show E = mean(skipmissing(mea[1].energy[N]));
@show δE = std(skipmissing(mea[1].energy[N]));

# =================================================================================================

# [det_C(0.2, 10, 1.0) for _ in 1:5]
# E = pmap(zeros(20), on_error=identity) do _
#     mea = main(M, T, N, C);
#     Eᵥ = mean(skipmissing(mea[1].energy_virial[N]))
#     Eₜ = mean(skipmissing(mea[1].energy[N]))
#     return [Eᵥ, Eₜ]
# end
# E = Matrix(reduce(hcat, E)')
# @show mean(E[:, 1])
# @show std(E[:, 1])
# # # energy(2, 1.0, 2)
# @show mean(E[:, 2])
# @show std(E[:, 2])

# =================================================================================================
# C = 4e-5
# τ = 0.05
# Ns = Int64[2, 3, 4, 5]
# Ts = 10 .^ [-1.2+0.2*k for k in 0:7]
# pmap(Iterators.product(Ns, Ts)) do (N, T)
#     mea = main(Int64(round(1/(τ*T))), T, N, C);
#     Eᵥ = collect(skipmissing(mea[1].energy_virial[N]))
#     Eₜ = collect(skipmissing(mea[1].energy[N]))
#     l = length(Eᵥ)
#     matrix = hcat(ones(l) .* N, ones(l) .* T, Eᵥ, Eₜ)
#     df =  DataFrame(matrix, ["N", "T", "energy_virial", "energy_thermo"])
#     CSV.write("trapped_bose_energy_open_N$(N)_T$(str(T)).csv", df)
#     return matrix
# end

# =================================================================================================

# Ms=Int64[10, 50, 100, 200, 300]
# Ts=Float64[5.0, 2.0, 1.0, 0.5, 0.1, 0.05, 0.01]
# E = pmap(Iterators.product(Ts, Ms), on_error=identity) do (T,M)
#     mea = main(M, T, N, C)
#     Eᵥ = mean(skipmissing(mea[1].energy_virial[N]))
#     Eₜ = mean(skipmissing(mea[1].energy[N]))
#     return [T, M, Eᵥ, Eₜ]
# end
# E = Matrix(reduce(hcat,E)')
# df = DataFrame(E, ["temperature", "timeslices", "energy_virial", "energy_thermo"])
# # insertcols!(df, 1, :temperature => Ts)
# # insertcols!(df, 2, :tau => τ .* ones(length(Ts)))
# @show df;
# CSV.write("trapped_bose_energy_N2_varying_temperature_and_timeslices.csv", df)

# =================================================================================================

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
# end.
