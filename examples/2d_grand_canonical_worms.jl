# potential = generate_V(cyclic_4, scale, depth; sym=:plot);
# x = LinRange(-L, L, 1000); y = LinRange(-L, L, 1000)
# X, Y = meshgrid(x, y); V = potential.(X, Y);
# plt = plot(legend = :none, colorbar=:none)
# surface!(X, Y, V, c=cgrad(:Greens_9), camera = (30, 20))


using Pimc
# using BenchmarkTools, OwnTime, StatsBase
using JLD2, FileIO, DelimitedFiles
using Plots
gr(show = true)
println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl");
include(pwd()*"/examples/tools/potentialtools.jl");
# include(pwd()*"/examples/tools/plottools.jl")
# include("tools\\savetools.jl")
println("Ended precompiling")
wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)
str(x::Float64) = replace(string(x), "."=>"_")

function main(L, g, M, T, n, times, potential)
    propint = build_prop_int(round(√(L^2+L^2), RoundUp), g, 1/(T*M))
    s = System(
        potential;
        dim = 2,
        M = M,
        λ = 1/(π^2),
        N = 5,
        L = L,
        T = T,
        C = 3e-5,
        g = g,
        μ = 0.3,
        worms = true,
        gc = true,
        propint = propint,
        interactions = true,
        length_measurement_cycle=1
        );
    info(s)


    Zupdates = [
        (1, SingleCenterOfMass(s, 1.0)),
        (1, ReshapeLinear(s, 5)),
        (1, ReshapeSwapLinear(s, 20))
    ]
    updates, Gupdates = add_worm_updates(s::System, Zupdates::Updates, m = 3)
    Gmeasurements = [
        NumbOfParticle(2*times*n*s.Ncycle)
    ]
    Zmeasurements = [
        SuperFluid(s, 2*times*n*s.Ncycle)
    ]

    # println("thermalization")
    # run!(s, times*n, updates)

    s.N_MC[0] = 0
    while s.N_MC[0] < n*times
        println(" ")
        @info s.N_MC[0]
        @info s.N
        @info get_C(s.C)
        @info acceptance(s.C.queue) # not N_g/N_z

        run!(s, n, Gupdates; Gmeasurements = Gmeasurements)
    end
    return Gmeasurements[1]
end
# propint = load("Store/Interactionpropagators/tau0_01/propint-L6_0g$(str(g))tau0_01.jld2", "prop_int");
# println("building interaction propagator")
cyclic_4 = [2π * k / 4 for k in 0:3];
scale = 1.0; depth = 5.0;
potential = generate_V(scale, depth, :l65)

n = 10_000; times = 10
L = 4.0
g = 16.16
M = 10
T = 1.0
mea = main(L, g, M, T, n, times, potential)
N = collect(skipmissing(mea.particle))
# mean(collect(skipmissing(mea[1].sf[5])))

plt = plot(fontfamily="Computer Modern",
    palette=:Set1_9,
    foreground_color_legend=nothing,
    background_color_legend=nothing,
    # legend_position=(0.17, 0.86),
    grid=false,
    legendfontsize=12,
    dpi=200,
    minorticks=10,
    size = (600, 200),
    framestyle = :box, legend = false
)
# plot!([10, 10])
plot!(N, linecolor=2)
savefig(plt, "number_of_particle.svg")





# @code_warntype main()
# @profview main()
# owntime(stackframe_filter=filecontains("Pimc.jl"))


# mutable struct HeadTimeSlice <: Measurement
#     headtimeslice::Vector{Union{Missing,Int64}}

#     HeadTimeSlice(n = 20_000) = new(Vector{Union{Missing,Int64}}(missing, n))
# end

# function (d::HeadTimeSlice)(s::System)
#     if s.worms > 0
#         ntail = rand(findall(i -> isa(i, Worm) && !isa(i.tail, Missing), s.world))
#         NpolW, polW = subcycle(s.world, ntail)
#         nhead = polW[NpolW]
#         d.headtimeslice[findfirst(ismissing, d.headtimeslice)] = s.world[nhead].head
#     end
# end
