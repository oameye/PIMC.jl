using Pimc
using StatsBase, LinearAlgebra
using Plots, LaTeXStrings
gr(show = true);
Plots.GRBackend();

println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
include(pwd()*"/examples/tools/plottools.jl")
# include("C:/Users/drandran/Documents/Pimc.jl/examples/tools/savetools.jl");
println("Ended precompiling")

function zero_worldline!(s::System)::Nothing
    r = zeros(Float64, s.M, s.dim)
    V = zeros(Float64, s.M)
    bins = zeros(Int64, s.M)
    for j in eachindex(V)
        bins[j] = bin(r[j, :], s.nbins, s.L)
        V[j] = lnV(r[j, :], r[mod1(j + 1, M), :], s.τ, s.V)
    end
    p = Particle(r, V, bins, s.N+1)
    push!(s.world, p)
    update_nnbins!(s)
    s.N += 1
    s.worms = 0
    nothing
end

function worldline!(s::System)::Nothing
    r = ones(Float64, s.M, s.dim) .* rand()
    V = zeros(Float64, s.M)
    bins = zeros(Int64, s.M)
    for j in eachindex(V)
        bins[j] = bin(r[j, :], s.nbins, s.L)
        V[j] = lnV(r[j, :], r[mod1(j + 1, M), :], s.τ, s.V)
    end
    p = Particle(r, V, bins, s.N+1)
    push!(s.world, p)
    update_nnbins!(s)
    s.N += 1
    s.worms = 0
    nothing
end

function init(M, N, T, C)::System
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
    for _ in 1:N
        worldline!(s)
    end
    # display(worldlines(s))
    return s
end

function linehistogram(E, plt; nbins=70, label=L"")
    h = normalize(fit(Histogram, E, nbins=nbins), mode=:density)
    r = h.edges[1]
    x = first(r)+step(r)/2:step(r):last(r)
    return plot!(plt, x, h.weights, label=label)
end

function filtermissing(v)
    Pimc.disallowmissing(filter(x->!ismissing(x), v))
end

M = 10
T = 1.0
C = 1.0
N = 2
n = 100_000
m = 5

s = init(M, N, T, C);

reshape_linear = ReshapeLinear(s, m);
reshape_swap = ReshapeSwapLinear(s, m);
reshape = [reshape_linear, reshape_swap];
Et_reshape = Energy(s, n);

for _ in 1:n
    reshape_linear(s)
end
pass = true
while pass
    rand(reshape)(s)
    Et_reshape(s)
    l = length(filtermissing(Et_reshape.energy[N]))
    if l >= n
        pass = false
    end
end

plt = plot(
    fontfamily="Computer Modern",
    palette=:Set1_9,
    grid=false,
    legend=true,
    legend_position=(0.42, 0.83))
plt = linehistogram(filtermissing(Et_reshape.energy[N]), plt, label="reshape")
display(plt)

s = init(M, N, T, C);

OC = Threshold(m, 1, s.M-1, 0.4, 0.6);
SW = NumbOfSlices(9, 1, s.M-1, 0.4, 0.6);
openc = OpenHeadC(OC);
closec = CloseHeadC(OC);
movehead = MoveHead(OC);
swap =  SwapC(SW);
moveswap = [movehead, swap];
Et_worm = Energy(s, n);

for _ in 1:n
    reshape_linear(s)
end
pass = true
while pass
    while s.worms == 0
        openc(s)
    end
    for _ in 1:5
       rand(moveswap)(s)
    end
    while s.worms > 0
        closec(s)
    end
    Et_worm(s)
    l = length(filtermissing(Et_worm.energy[N]))
    if l >= n
        pass = false
    end
end
length(filtermissing(Et_worm.energy[N]))
plt = linehistogram(Pimc.disallowmissing(Et_worm.energy[N]), plt, label="worms")
plot(plt, legend_position=(0.22, 0.93), xlims=(-20, 20))

# s = init(M, T, C);

# opengc = OpenGC(OC);
# closegc = CloseGC(OC);
# Ev_oc_gc = VirialEnergy(s, n);
# Et_oc_gc = Energy(s, n);

# for _ in 1:n
#     reshape_linear(s)
# end
# for _ in 1:n
#     pass = true
#     while pass
#         for _ in 1:5
#             opengc(s)
#             closegc(s)
#         end
#         pass = s.worms > 0 ? true : false
#     end
#     # Ev_reshape(s)
#     Et_oc_gc(s)
# end

# plt = linehistogram(Pimc.disallowmissing(Et_oc_gc.energy[N]), plt, label="open-close grand canonical")
# plot(plt, legend_position=(0.22, 0.93), xlims=(-200, 50))
