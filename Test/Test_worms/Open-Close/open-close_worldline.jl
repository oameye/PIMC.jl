using Pimc
using StatsBase
using Plots, LaTeXStrings
gr(show = true);
Plots.GRBackend();

println("Compiling plot and save function")
include("C:/Users/drandran/Documents/Pimc.jl/examples/tools/logtools.jl");
include("C:/Users/drandran/Documents/Pimc.jl/examples/tools/plottools.jl");
# include("C:/Users/drandran/Documents/Pimc.jl/examples/tools/savetools.jl");
println("Ended precompiling")

function zero_worldline!(s::System)::Nothing
    r = zeros(Float64, s.M, s.dim)
    V = zeros(Float64, s.M)
    bins = zeros(Float64, s.M)
    p = Particle(r, V, bins, 1)
    push!(s.world, p)
    update_nnbins!(s)
    s.N = 1
    s.worms = 0
    nothing
end

function init(M, T, C)::System
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
                N = 0,
                L = 100.0,
                T = T,
                SC = SC,
                worms = true,
                gc = false,
                modifyC = false,
                measure_scheme = :c,
                length_measurement_cycle = 1
        );

    zero_worldline!(s)
    # display(worldlines(s))
    return s
end

M = 100
T = 1.0
C = 1.0
n = 10_000
m = 20
s = init(M, T, C);
function nextframe(s, m, n)
    reshape_linear = ReshapeLinear(s, m);
    for _ in 1:n
        reshape_linear(s)
    end
    plt = worldlines(s, linecolor = :black)

    OC = Threshold(m, 1, s.M-1, 0.4, 0.6)
    openc = OpenHeadC(OC);
    closec = CloseHeadC(OC);
    for _ in 1:n
        openc(s)
        closec(s)
    end
    while s.worms > 0
        closec(s)
    end
    plt = worldlines(s, plt, linecolor = :blue)
end
display(nextframe(s, m, n))

