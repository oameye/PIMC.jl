using Pimc
# using BenchmarkTools
# using JLD2, FileIO
using Plots
# using StatsBase
# using OwnTime
gr(show = true)
Plots.GRBackend()

println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
include(pwd()*"/examples/tools/plottools.jl")
# include("tools\\savetools.jl")
println("Ended precompiling")

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

s = System(
    (r) -> 0.0;
    dim = 2,
    M = 100,
    N = 2,
    L = 4.0,
    T = 1.0);

n = 100_000
Zupdates = [
    (1, PolymerCenterOfMass(s, 1.0)),
    (1, ReshapeLinear(s, 20))
    # (1, ReshapeSwapLinear(s, 3))
];
# updates, Gupdates = add_worm_updates(s, Zupdates);

info(s)
# display(worldlines(s))
# println("running")
# function nextframe(s, updates)
#     run!(s, 1, updates)
#     return worldlines(s)
# end
# anim = @animate for _ in 1:n
#     nextframe(s, updates)
# end
# mp4(anim, "test.mp4", fps = 50)
run!(s, 10_000, Zupdates)
for _ in 1:n
    # wait_for_key("Press enter to continue")
    run!(s, 1, Zupdates)

    display(worldlines(s; dim = 1))
    wait_for_key("Press enter to continue")
end
# end
# main()

# @code_warntype main()
# @profview main()
# owntime(stackframe_filter=filecontains("Pimc.jl"))







