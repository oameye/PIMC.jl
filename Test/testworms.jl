using StatsBase

# if 2*s.L/s.nbins < s.a
#     @info "binsize ($(2*s.L/s.nbins)) is smaller than 2d scattering length ($(s.a))"
#     return
# end

function runworms!(s::System, n::Int64, updates::Updates)::Nothing
    (every, fs) = Pimc.unzip(updates)
    weights = Weights(1 ./ every)
    for _ = 1:n
        Pimc.apply!(s, fs, weights)
        permerror(s)
        checkworms(s)
    end
    nothing
end

function debugpermlist(s::System)
    perm = Vector{String}(undef, length(s.world))
    for (i,p) in enumerate(s.world)
        if isa(p, Pimc.Worm)
            if !isa(p.tail, Missing) && !isa(p.head, Missing)
                perm[i] = "W"
            elseif isa(p.tail, Missing)
                perm[i] = "H"
            elseif isa(p.head, Missing)
                perm[i] = "T$(p.next)"
            else
                @error "Not a worm"
            end
        else
            perm[i] = "$(p.next)"
        end
    end
    return perm
end

@testset "Worms" begin


scale = 0.5
depth = 8.

v1d(x::Vector{Float64}) = depth * sin(2π * x[1] * scale)^2

s = System(v1d; dim=1, M=100, N=5, L=4., T= 1., μ=3.55, worms=true)

worms = [
    (1, Insert(s, 25)),
    (1, Remove(s, 100)),
    (1, Open(s, 25)),
    (1, Close(s, 100)),
    (1, Advance(s, 40)),
    (1, Recede(s, 40)),
    (1, Swap(s, 5))
]

# updates = s.worms ? append!(updates, worms) : updates

n = 1
times = 10_000

for i in 1:times
    runworms!(s, n, worms)

    perm = debugpermlist(s)
    for l in eachindex(s.world)
        if count(==("$(l)"), perm) > 1
            @test error("permutations are not consistant\n $(l) is more than one in $(perm)")
            break
        end
        if count(==("H"), perm) > 1
            @test error("permutations are not consistant\n H is more than one in $(perm)")
            break
        end
        if count(==("W"), perm) > 1
            @test error("permutations are not consistant\n W is more than one in $(perm)")
            break
        end
        if !occursin("H", perm[l]) && !occursin("W", perm[l])
            if occursin("T", perm[l])
                if parse(Int64, perm[l][2:end]) > length(s.world)
                    @test error("permutations are not consistant\n $(parse(Int64, perm[l])) at index $(l) is out of bound in $(perm)")
                end
            else
                if parse(Int64, perm[l]) > length(s.world)
                    @test error("permutations are not consistant\n $(parse(Int64, perm[l])) at index $(l) is out of bound in $(perm)")
                end
            end
        end
        N = count(i -> !occursin("H", i) && !occursin("T", i)  && !occursin("W", i), perm)
        if N != s.N
            @test error("The number of particles in stored in the struct system is $(s.N) but I count $(N)")
        end
        a = count(i -> isa(i, Pimc.Worm), s.world)
        if s.gsec != a
            @test error("The number of Worms in stored in the struct system is $(s.gsec) but I count $(worms)")
        end

        if !isempty(findall(i -> isa(i, Pimc.Worm) && !isa(i.tail, Missing), s.world))
            ntail = rand(findall(i -> isa(i, Pimc.Worm) && !isa(i.tail, Missing), s.world))
            NpolW, polW = subcycle(s.world, ntail)
            nhead = polW[NpolW]
            if ismissing(s.world[nhead].head)
                @test error("The time slice head of is missing")
            elseif ismissing(s.world[ntail].tail)
                @test error("The time slice tail of is missing")
            end
            if ismissing(s.world[nhead].r[s.world[nhead].head, 1])
                @test error("The time slice head in the position vector is missing")
            elseif ismissing(s.world[ntail].r[s.world[ntail].tail, 1])
                @test error("The time slice tail in the position vector is missing")
            end
            if ntail==nhead && s.gsec !=1
                @test error("only one worm but gsec = $(s.gsec)")
            end
            if ntail!=nhead && s.gsec !=2
                @test error("two worms but gsec = $(s.gsec)")
            end
        end
    end
end

end # @testset "Worms"