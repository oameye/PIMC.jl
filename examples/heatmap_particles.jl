using Distributed
addprocs(6)

@everywhere begin
using Pimc
using JLD2, FileIO
using StatsBase

println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
include(pwd()*"/examples/tools/plottools.jl")
# include("tools\\savetools.jl")
println("Ended precompiling")

# wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)
str(x::Float64) = replace(string(x), "."=>"_")

function newsystem(g::Float64, μ::Float64)::System
    propint = load("Store\\Interactionpropagators\\tau0_01\\propint-L6_0g$(str(g))tau0_01.jld2", "prop_int")
    System(
        _ -> 0.0;
        dim = 2,
        M = 100,
        N = 5,
        L = 4.0,
        T = 1.0,
        g = g,
        μ = μ,
        worms = true,
        propint = propint,
        interactions = true)
end
function main(g, μ)
    # @show (g, μ)
    s = newsystem(g, μ)
    worms = worm_updates(s)
    n = 1_000_000
    measurements = [
        (1, NumbOfParticle(n))
    ]
    run!(s, n, worms, measurements)
    particles = mean(collect(skipmissing(measurements[1][2].particle)))
    # @show (g, μ, particles)
    return (g, μ, particles)
end

end
gs = 0.16:1.0:16.16
μs = range(2.0, 4.0, 16)
status = pmap(Iterators.product(gs, μs), on_error=identity) do (g, μ)
    main(g, μ)
end

open("heatmap_particles.txt", "w") do file
    for out in status
        @show status
        println(file, out)
    end
end
