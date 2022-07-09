using CSV
using DataFrames
#TODO worms compatiable saving
function paths_df(s::System)
    if s.dim == 2
        names = collect(Iterators.flatten(zip("p" .* string.(1:s.N) .* " x", "p" .* string.(1:s.N) .* " y")))
    else
        names = "p" .* string.(1:s.N) .* " x"
    end
    pushfirst!(names, "tau")
    data = [zeros(s.M + 1) for _ in 1:(s.dim*s.N+1)]
    data[1] = [j * s.β / s.M for j in 0:s.M]
    if s.dim == 2
        dummy = -1
    else
        dummy = 0
    end
    for n in 1:s.N
        Npol, pol = subcycle(s.world, n)
        if s.dim == 2
            dummy += 1
        end
        for d in 1:s.dim
            data[n+d+dummy] = push!(s.world[n].r[:, d], s.world[pcycle(s.M + 1, pol, Npol, s.M)].r[1, d])
        end
    end
    return DataFrame(data, names)
end

function save_paths(s::System)
    df = paths_df(s)
    CSV.write("paths.csv", df)
    println("saved paths")
end

function density_df(s::System, d::Density)
    den = d.dens ./ (d.bin * d.ndata)
    if den isa Vector
        df = DataFrame(:Density => den)
    else
        df = DataFrame(den, :auto)
    end
    pos = range(-s.L, s.L, length = d.nbins)
    insertcols!(df, 1, :pos => pos)
    return df
end

# function save_density(s::System, d::Density)
#     df = density_df(s, d)
#     CSV.write("densityM$(s.M)beta$(s.β)N$(s.N).csv", df)

#     # sanity check
#     numb = (sum(d.dens) / d.ndata)
#     if abs(numb - s.N) > 1e-4
#         println("Error save_density: integrated density = $numb, should be $(s.N) if x-range is ok")
#         # exit(0)
#     end
#     println("saved density")
# end
function save_density(s::System, d::Density, g, name, V)
    df = density_df(s, d)
    CSV.write("density$(name)M$(s.M)beta$(s.β)N$(s.N)g$(g)V$(V).csv", df)

    # sanity check
    # numb = (sum(d.dens) / d.ndata)
    # if abs(numb - s.N) > 1e-4
    #     println("Error save_density: integrated density = $numb, should be $(s.N) if x-range is ok")
    #     # exit(0)
    # end
    println("saved density")
end