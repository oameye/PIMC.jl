using Pimc
using JLD2, FileIO


println("Compiling plot and save function")
include("tools/logtools.jl")
include("tools/plottools.jl")
include("tools/savetools.jl")
println("Ended precompiling")

Coord = Vector{Float64}
Endpoint = Vector{Coord}

function create_config!(s::System, config::Vector{Endpoint}, perm::Vector{Int64})
    ps = Vector(undef, length(config))
    for n in eachindex(ps)
        V = Vector{Float64}(undef, s.M)
        r = Matrix{Float64}(undef, s.M, s.dim)
        r[begin, :] = config[n][1]
        r[end, :] = config[n][2]
        levy!(r, s.τ, s.L)

        pass = n == 1 ? false : true
        while pass
            @label escape_label
            for m in 1:s.M
                for i in 1:(n-1)
                    d = norm(distance.(ps[i].r[m, :], r[m, :], s.L))
                    if d < s.a
                        levy!(r, s.τ, s.L)
                        @goto escape_label
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
    s.world = ps
end

function test_reshape_config(s::System)
    config = [[Coord([3.5, 0]), Coord([-3.2, 0])]]
    perm = Int64[1]
    create_config!(s, config, perm)
    update_nnbins!(s)
    s.N = 1
end

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

function main()
    propint = load("Store\\Interactionpropagators\\tau0_01\\propint-L6_0g1_16tau0_01.jld2", "prop_int")

    s = System(_ -> 0; M = 30, N = 1, L = 4.0, T = 1, g = 1.16, propint = propint, interactions = false)

    info(s) # Layout with general information for the initiation of the simulation

    acc = false
    while !acc
        test_reshape_config(s)
        u = s
        display(worldlines(s))
        # savefig("Store\\fig\\plot$i")
        wait_for_key("press any key to continue")
        acc = ReshapeLinear(s, 100)(s)
        display(worldlines(s))
        println("debug")
    end

end
println("Ended precompiling")
for _ in 1:1
    main()
end