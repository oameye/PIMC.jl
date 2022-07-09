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

Coord = Vector{Union{Float64,Int64}}
Endpoint = Vector{Coord}

function create_config!(s::System, config::Vector{Endpoint}, perm::Vector{Int64}, worm_config::Vector{Endpoint})
    ps = Vector(undef, length(config) + length(worm_config))
    for n in eachindex(config)
        V = Vector{Float64}(undef, s.M)
        r = Matrix{Float64}(undef, s.M, s.dim)
        r[begin, :] = config[n][1]
        r[end, :] = config[n][2]
        levy!(r, s.τ, s.L)

        pass = n == 1 ? false : true
        while pass
            @label escape_label_particle
            for m in 1:s.M
                for i in 1:(n-1)
                    d = norm(distance.(ps[i].r[m, :], r[m, :], s.L))
                    if d < s.a
                        levy!(r, s.τ, s.L)
                        @goto escape_label_particle
                    end
                end
            end
            pass = false
        end

        # Pre-calc propagator for each time slice
        for j in eachindex(V)
            V[j] = s.lnV(r[j, :], r[mod1(j + 1, s.M), :])
        end
        ps[n] = Particle(r, V, zeros(s.M), perm[n])
    end

    #adds worms from the ground up, i.e. 1 to M
    V = Vector{Union{Float64,Missing}}(missing, s.M)
    r = Matrix{Union{Float64,Missing}}(missing, s.M, s.dim)
    for n in eachindex(worm_config)
        m = Int64(worm_config[n][2][3] - worm_config[n][1][3])
        r′ = zeros(Float64, m + 1, s.dim)
        r′[begin, :] = worm_config[n][1][1:s.dim]
        r′[end, :] = worm_config[n][2][1:s.dim]
        levy!(r′, s.τ, s.L)

        pass = n == 1 ? false : true
        while pass
            @label escape_label_worms
            for m in eachindex(r′[:, 1])
                for i in eachindex(config)
                    d = norm(distance.(ps[i].r[Int64(worm_config[n][1][3] + m - 1), :], r′[m, :], s.L))
                    if d < s.a
                        levy!(r′, s.τ, s.L)
                        @goto escape_label_worms
                    end
                end
            end
            pass = false
        end
        r[Int64(worm_config[n][1][3]):Int64(worm_config[n][2][3]), :] = r′

        # Pre-calc propagator for each time slice
        for j in eachindex(r′[:, 1])
            m = Int64(worm_config[n][1][3] + j - 1)
            V[m] = 0.0
        end

        if n == 1 && length(worm_config) > 1
            ps[length(config)+n+1] = Worm(r, V, zeros(s.M), missing, Int64(worm_config[n][2][3]), missing)
        elseif n == 2
            ps[length(config)+n-1] = Worm(r, V, zeros(s.M), Int64(worm_config[n][1][3]), missing, length(config) + n)
        else
            ps[length(config)+n] = Worm(r, V, zeros(s.M), Int64(worm_config[n][1][3]), Int64(worm_config[n][2][3]), missing)
        end
    end
    s.world = ps
end

function test_advance_config(s::System)
    config = Vector{Endpoint}()
    perm = Int64[]
    worm_config = [[Coord([-1.1, -1.1, Int64(30)]), Coord([-1.1, -1.1, Int64(45)])]]
    create_config!(s, config, perm, worm_config)
    update_nnbins!(s)
    s.N = 0
    s.worms = length(worm_config)
    s.worm_algorithm = true
end

function main()
    s = System(_ -> 0;
        dim = 2,
        M = 100,
        N = 2,
        L = 4.0,
        T = 1.0,
        gc = false,
        worms = true)

    AR = Threshold(10, 1, s.M-1, 0.4, 0.6)

    for _ in 1:1
        test_advance_config(s)
        @info "new config done"
        info(s)
        display(worldlines(s))
        acc = MoveTail(AR)(s)
        @info "MoveTail accepted: $acc"
        display(worldlines(s))
    end
    # @info remove_accepted/100
end
main()