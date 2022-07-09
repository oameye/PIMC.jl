# write nn with CartesianIndices and SMatrix

"""
	bin_index(x::Int64, y::Int64, w::Int64)::Int64

Return the index of the bin in the `nn`` vector given the indices `x` and `y` of the 2D spatial grid. The argument `w` is the number of bins of the 2D spatial grid in one dimension.
"""
function bin_index(x::Int64, y::Int64, w::Int64)::Int64
    x + w * y + 1
end
"""
	bin_index(x::Int64 w::Int64)

Return the index of the bin in the `nn` vector given the indices `x` of the 2D spatial grid. The argument `w` is the number of bins of the spatial grid in one dimension.
"""
function bin_index(x::Int64, w::Int64)::Int64
    x + 1
end

"""
	bin(r::Vector{Float64}, nbins::Int64, L::Float64)::Int64

Return the index of the bin in the `nn` vector given the particle position `r`, the number of bins of the spatial grid in one dimension `nbins`, and the width of the system(lattice) `2*L`.
"""
# Extract the whole part of the particle coordinates to assign a grid position
function bin(r::Vector{Float64}, nbins::Int64, L::Float64)::Int64
    ibin = floor.(Int64, (r .+ L) ./ (2 * L / nbins))
    bin_index(ibin..., nbins)
end
function bin(r::SubArray{Float64,1}, nbins::Int64, L::Float64)::Int64
    ibin = floor.(Int64, (r .+ L) ./ (2 * L / nbins))
    bin_index(ibin..., nbins)
end

"""
	psort!(ps::Vector{Worldline}, bins::Vector{Vector{Int64}}, nbins::Int64, j::Int64, L::Float64)::Nothing

For a give timeslice `j` the function goes over every particle in the `world`, computes which bin it the particle finds himself and add the bin to the bins vector `bins`. The bin is also added the bin vector in the `Worldline` struct of `world`.
"""
function psort!(world::Vector{Worldline}, bins::Vector{Vector{Int64}}, nbins::Int64, j::Int64, L::Float64)::Nothing
    for (i, p) in enumerate(world)
        if !ismissing(p.r[j, 1])

            b = bin(disallowmissing(p.r[j, :]), nbins, L)
            push!(bins[b], i)

            # Also link bin index to particle
            p.bins[j] = b
        end
    end
	nothing
end

# Find adjacent cells given periodic boundary conditions
function bin_neighbors(b::Int64, nbins::Int64, dim::Int64)::Vector{Int64}
    if dim == 2
        x = rem(b - 1, nbins)
        y = div(b - 1, nbins)
        [bin_index(mod(x + dx, nbins), mod(y + dy, nbins), nbins)
         for (dx, dy) in ((0, 0), (-1, 1), (0, 1), (1, 1), (-1, 0), (1, 0), (-1, -1), (0, -1), (1, -1))]
    else
        x = rem(b - 1, nbins)
        [bin_index(mod(x + dx, nbins), nbins) for dx in (0, -1, 1)]
    end
end

"""
	within_neigbourhood(s::System, i::Int64, j::Int64; exception = [])

Finds and collect the particles in the bin and the neigbouring bins of particle `i` at time slice `j`. It expand the neigbouring bins range until the output is not empty and excludes the beads at timeslice `j` of the particles in the list `exceptions`.
"""
function within_neigbourhood(s::System, i::Int64, j::Int64; exceptions=Int64[])::Vector{Int64}
    if ismissing(s.world[i].bins[j])
        error("bin in particle struct $i at time slice $j is missing\n
        whereas the position $(s.world[i].r[j, :])")
    end
    p_nnbins = filter!(x -> x ∉ exceptions, vcat(s.nn[j][s.nbs[s.world[i].bins[j]]]...))
    return p_nnbins
end

"""
	within_neigbourhood(s::System, i::Int64, j::Int64; exception = [])

Finds and collect the particles in the bin and the neigbouring bins of position `r` at time slice `j`. It expand the neigbouring bins range until the output is not empty and excludes the beads at timeslice `j` of the particles in the list `exceptions`.
"""
function within_neigbourhood(s::System, r::Vector{Float64}, j::Int64; exceptions=Int64[])::Vector{Int64}
    b = bin(r, s.nbins, s.L)
    if b < 0
        @show r
        @show b
    end
    p_nnbins = filter!(x -> x ∉ exceptions, vcat(s.nn[j][s.nbs[bin(r, s.nbins, s.L)]]...))
    return p_nnbins
end

"""
	within_neigbourhood(nn::Vector{Vector{Vector{Int64}}}, nbs::Vector{Vector{Int64}}, b::Int64, j::Int64)

Finds and collect the particles in the bin `b` and the neigbouring bins `nbs`` at time slice `j` all stored in `nn`. It expand the neigbouring bins range until the output is not empty and excludes the beads at timeslice `j` of the particles in the list `exceptions`.
"""
function within_neigbourhood(nn::Vector{Vector{Vector{Int64}}}, nbs::Vector{Vector{Int64}}, b::Int64, j::Int64; exceptions=Int64[])::Vector{Int64}
    p_nnbins = filter!(x -> x ∉ exceptions, vcat(nn[j][nbs[b]]...))
    return p_nnbins
end

function find_nns(s::System, i::Int64, j::Int64; exceptions=Int64[], leafsize=1)::Vector{Int64}
    # collect the particles within cell and it neigbours at a certain timeslice j
    p_nnbins = within_neigbourhood(s, i, j, exceptions=exceptions)
    if isempty(p_nnbins)
        return p_nnbins
    end
    # Combine coordinates to work with KDTree
    # FIXME error disallowmissing (very rare)
    for n in p_nnbins
        if ismissing(s.world[n].r[j, 1])
            error("position of particle $n at timeslice $j is missing
            \n the vector bin gives $(s.world[n].bins[j])
            \n and nn grid bin contains the particles $(filter!(x -> x ∉ exceptions, vcat(s.nn[j][s.nbs[s.world[i].bins[j]]]...)))
            ")
        end
    end
    coords_particles = hcat([disallowmissing(s.world[n].r[j, :]) for n in p_nnbins]...) .+ s.L
    target_coord = s.world[i].r[j, :] .+ s.L

    # Make tree and get NNs
    balltree = BallTree(coords_particles, PeriodicEuclidean(zeros(s.dim) .+ (2 * s.L)), leafsize=leafsize)
    idxs = inrange(balltree, target_coord, (2 * s.L) / s.nbins)
    return [p_nnbins[n] for n in idxs]
end

function find_nns(s::System, r::Vector{Float64}, j::Int64; exceptions=Int64[], leafsize=1)::Vector{Int64}
    # collect the particles within cell and it neigbours at a certain timeslice m
    p_nnbins = within_neigbourhood(s, r, j, exceptions=exceptions)
    if isempty(p_nnbins)
        return p_nnbins
    end
    # Combine coordinates to work with KDTree
    # FIXME error disallowmissing (very rare)
    for n in p_nnbins
        if ismissing(s.world[n].r[j, 1])
            error("position of particle $n at timeslice $j is missing
            \n the vector bin gives $(s.world[n].bins[j])
            \n and nn grid bin contains the particles $(filter!(x -> x ∉ exceptions, vcat(s.nn[j][s.nbs[bin(r, s.nbins, s.L)]]...)))
            ")
        end
    end
    coords_particles = disallowmissing(hcat([s.world[n].r[j, :] for n in p_nnbins]...) .+ s.L)
    target_coord = r .+ s.L

    # Make tree and get NNs
    balltree = BallTree(coords_particles, PeriodicEuclidean(zeros(s.dim) .+ (2 * s.L)), leafsize=leafsize)
    idxs = inrange(balltree, target_coord, (2 * s.L) / s.nbins)
    return [p_nnbins[n] for n in idxs]
end

function find_nn(s::System, r::Vector{Float64}, j::Int64; exceptions = Int64[], leafsize = 1)::Int64
    # Select the particles within cell and it neigbours at a certain timeslice m
    p_nnbins = within_neigbourhood(s, r, j, exceptions=exceptions)
    if isempty(p_nnbins)
        return -1
    end

    # Combine coordinates to work with KDTree
    for n in p_nnbins
        if ismissing(s.world[n].r[j, 1])
            error("position of particle $n at timeslice $j is missing
            \n the vector bin gives $(s.world[n].bins[j])
            \n and nn grid bin contains the particles $(filter!(x -> x ∉ exceptions, vcat(s.nn[j][s.nbs[bin(r, s.nbins, s.L)]]...)))
            ")
        end
    end
    coords_particles = disallowmissing(hcat([s.world[n].r[j, :] for n in p_nnbins]...) .+ s.L)
    target_coord = r .+ s.L

    # Make tree and get NNs
    balltree = BallTree(coords_particles, PeriodicEuclidean(zeros(s.dim) .+ (2 * s.L)), leafsize=leafsize)
    n, _ = nn(balltree, target_coord)
    return p_nnbins[n]
end

# instead of making new structre just go over all beads and delete or replace
function update_nnbins!(s::System)
    for m in eachindex(s.nn)
        for b in eachindex(s.nn[m])
            s.nn[m][b] = Int64[]
        end
        for (n, p) in enumerate(s.world)
            if !ismissing(p.r[m, 1])
                b′ = bin(disallowmissing(p.r[m, :]), s.nbins, s.L)
                push!(s.nn[m][b′], n)

                p.bins[m] = b′
            end
        end
    end
end

function update_nn_bead!(s::System, n::Int64, j::Int64, r′::Vector{Float64})
    filter!(x -> x != n, s.nn[j][s.world[n].bins[j]])
    b = bin(r′, s.nbins, s.L)
    push!(s.nn[j][b], n)
    s.world[n].bins[j] = b
end
function update_nn_bead!(s::System, n::Int64, j::Int64, r′::SubArray{Float64,1})
    filter!(x -> x != n, s.nn[j][s.world[n].bins[j]])
    b = bin(r′, s.nbins, s.L)
    push!(s.nn[j][b], n)
    s.world[n].bins[j] = b
end

function add_nn!(nn::NearestNeighbors, pol::Vector{Int64}, w::Vector{Worldline})
    for n in pol
        for (j, b) in pairs(skipmissing(w[n].bins))
            push!(nn[j][b], n)
        end
    end
end
function add_nn!(nn::NearestNeighbors, n::Int64, bins::Vector{Union{Int64,Missing}})
    for (j, b) in pairs(skipmissing(bins))
        push!(nn[j][b], n)
    end
end
function add_nn!(nn::NearestNeighbors, n::Int64, bins::Vector{Int64})
    for (j, b) in pairs(bins)
        push!(nn[j][b], n)
    end
end

function rm_nn!(nn::NearestNeighbors, pol::Vector{Int64}, w::Vector{Worldline})
    for n in pol
        for (j, b) in pairs(skipmissing(w[n].bins))
            filter!(x -> x != n, nn[j][b])
        end
    end
end

function rm_nn!(nn::NearestNeighbors, n::Int64, bins::Vector{Union{Int64,Missing}})
    for (j, b) in pairs(skipmissing(bins))
        filter!(x -> x != n, nn[j][b])
    end
end
function rm_nn!(nn::NearestNeighbors, n::Int64, bins::Vector{Int64})
    for (j, b) in pairs(bins)
        filter!(x -> x != n, nn[j][b])
    end
end

function fixindexing_nn!(nn::NearestNeighbors, n::Int64, n′::Int64, bins::Vector{Union{Int64,Missing}})
    for (j, b) in pairs(skipmissing(bins))
        filter!(x -> x != n, nn[j][b])
        push!(nn[j][b], n′)
    end
end
function fixindexing_nn!(nn::NearestNeighbors, n::Int64, n′::Int64, bins::Vector{Int64})
    for (j, b) in pairs(skipmissing(bins))
        filter!(x -> x != n, nn[j][b])
        push!(nn[j][b], n′)
    end
end



