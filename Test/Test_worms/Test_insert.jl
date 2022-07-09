using Pimc
using LinearAlgebra
# using BenchmarkTools
using JLD2, FileIO
using Plots
using OwnTime
gr(show = true)
Plots.GRBackend()

println("Compiling plot and save function")
include("c:\\Users\\drandran\\Documents\\Pimc.jl\\examples\\logtools.jl")
include("c:\\Users\\drandran\\Documents\\Pimc.jl\\examples\\plottools.jl")
include("c:\\Users\\drandran\\Documents\\Pimc.jl\\examples\\savetools.jl")
println("Ended precompiling")


function main()
    potential = x::Vector{Float64} -> 0

    propint = load("Interactionpropagators\\tau0_01\\propint-L6_0g1_16tau0_01.jld2", "prop_int")
    s = System(potential; dim = 2, M = 100, N = 5, L = 4.0, T = 1, g=1.16, μ = 2.0, worms = true, propint = propint)

    IR = Threshold(50, 0, 0, 1, s.M, 0.4, 0.6)

    config = [[[-1.0, -1.0], [1., 1.]], [[1., 1.], [-1., -1.]]]
    perm = [2, 1]
    create_config!(s, config)
    info(s)
    display(worldlines(s))
    Insert(IR)(s)
    display(worldlines(s))
end
main()