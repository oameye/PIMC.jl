using Distributed
addprocs(3)

using Plots, LaTeXStrings
gr(show = true);
Plots.GRBackend();
function linehistogram(E, plt; nbins=100, label)
    h = normalize(fit(Histogram, E, nbins=nbins))
    r = h.edges[1]
    x = first(r)+step(r)/2:step(r):last(r)
    return plot!(plt, x, h.weights, label=label)
end

using CSV, DataFrames

@everywhere begin
using Pimc
using StatsBase, LinearAlgebra

println("Compiling plot and save function")
include(pwd()*"/examples/tools/logtools.jl")
# include(pwd()*"/examples/tools/plottools.jl")
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

function filtermissing(v)
    Pimc.disallowmissing(filter(x->!ismissing(x), v))
end
function filtermissing!(v)
    Pimc.disallowmissing(filter!(x->!ismissing(x), v))
end

const M = 20
const T = 1.0
const C = 1e-2
const N = 2
const n = 1_000_000
const m = 5

function Ereshape(M, T, C, N, n, m)
    s = init(M, N, T, C);
    reshape_linear = ReshapeLinear(s, m);
    # reshape_swap = ReshapeSwapLinear(s, m);
    # reshape = [reshape_linear, reshape_swap];
    Et_reshape = Energy(s, n*2);
    @info "thermalize reshape"
    for _ in 1:n
        reshape_linear(s)
    end
    @info "reshape"
    global pass1 = true
    while pass1
        reshape_linear(s)
        Et_reshape(s)
        l = length(filtermissing(Et_reshape.energy[N]))
        if l >= n
            global pass1 = false
        end
    end
    return filtermissing(Et_reshape.energy[N])
end

function Eoc(M, T, C, N, n, m)
    s = init(M, N, T, C);
    reshape_linear = ReshapeLinear(s, m);
    OC = Threshold(m, 1, s.M-1, 0.4, 0.6);
    AR = Threshold(3, 1, s.M-1, 0.4, 0.6);
    openc = OpenGC(OC);
    closec = CloseGC(OC);
    oc = [openc, closec]
    Et_oc = Energy(s, n*2);
    @info "thermalize"
    for _ in 1:n
        reshape_linear(s)
    end
    @info "open and close"
    global pass1 = true
    while pass1
        for i in 1:5
            rand(oc)(s)
        end
        while s.worms > 0
            rand(oc)(s)
        end
        Et_oc(s)
        l = length(filtermissing(Et_oc.energy[N]))
        if l >= n
            global pass1 = false
        end
    end
    E_oc = filtermissing(Et_oc.energy[N])
    filter!(x -> x > -21, E_oc)
    return E_oc
end


function Ear(M, T, C, N, n, m)
    s = init(M, N, T, C);
    reshape_linear = ReshapeLinear(s, m);
    OC = Threshold(m, 1, s.M-1, 0.4, 0.6);
    AR = Threshold(3, 1, s.M-1, 0.4, 0.6);
    openc = OpenGC(OC);
    closec = CloseGC(OC);
    oc = [openc, closec]
    advance = Advance(AR);
    recede = Recede(AR);
    ar = [advance, recede];
    Et_worm = Energy(s, n*2);
    @info "thermalize"
    for _ in 1:n
        reshape_linear(s)
    end
    @info "advance and recede"
    global pass2 = true
    while pass2
        while s.worms == 0
            openc(s)
        end
        for _ in 1:5
           rand(ar)(s)
        end
        while s.worms > 0
            rand(ar)(s)
            rand(oc)(s)
        end
        Et_worm(s)
        l = length(filtermissing(Et_worm.energy[N]))
        if l >= n
            global pass2 = false
        end
    end
    E_ar = filtermissing(Et_worm.energy[N])
    length(filtermissing(Et_worm.energy[N]))
    filter!(x -> x > -21, E_ar)
    return E_ar
end

end #everywhere
E_r = remotecall(Ereshape, 2, M, T, C, N, n, m)
E_oc = remotecall(Eoc, 3, M, T, C, N, n, m)
E_ar = remotecall(Ear, 4, M, T, C, N, n, m)

E_r = fetch(E_r)
E_oc = fetch(E_oc)
E_ar = fetch(E_ar)

CSV.write("E_r.csv", DataFrame(E_r=E_r))
CSV.write("E_oc.csv", DataFrame(E_oc=E_oc))
CSV.write("E_ar.csv", DataFrame(E_ar=E_ar))

plt = plot(
    fontfamily="Computer Modern",
    palette=:Set1_9,
    grid=false,
    legend=true,
    legend_position=(0.42, 0.83),
    # yticks=false,
    xlabel="thermodynamic energy estimator",
    ylabel="normalized frequency",
    xlims=(-20, 20)
    )
plt = linehistogram(E_r, plt, label="reshape", nbins=70)
plt = linehistogram(E_oc, plt, label="OC", nbins=70)
plt = linehistogram(E_ar, plt, label="OC+AR*", nbins=70)
# plot!(plt, legend_position=(0.22, 0.93), xlims=(-20, 20))
savefig(plt, "energy_histogram_better.svg")

