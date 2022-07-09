# Real time plotting
using Plots, LaTeXStrings
gr(show = true)

Plots.GRBackend()
using CSV
using DataFrames
using Dates


# Find out if a worldline portals past into boundary at a timestep
function portal(wl::Vector{Float64}, boundary::Float64)::Vector{Int64}
    [i for i in 2:length(wl) if abs(wl[i] - wl[i-1]) > boundary]
end

function portal(wl::Vector{Union{Missing,Float64}}, boundary::Float64)::Vector{Int64}
    indices = []
    for i in 2:length(wl)
        if !ismissing(wl[i]) && !ismissing(wl[i-1])
            if abs(wl[i] - wl[i-1]) > boundary
                push!(indices, i)
            end
        end
    end
    indices
end

function pairwise(a::Vector{Int64}; weight = 1)::Vector{Vector{Int64}}
    pairs = []
    for (i, j) in zip(a[1:end-1], a[2:end])
        if (j - i) > weight
            push!(pairs, [i, j - 1])
        end
    end
    pairs
end

# Extract particle worldlines to xyz data
function extract_worldlines(s::System)
    plot_data = []
    t = collect(1:s.M)
    for p in s.world
        idx = []
        if s.dim == 2
            x, y = p.r[:, 1], p.r[:, 2]
            # Get all indices where a worldline crosses the box boundaries
            idx = [[1], portal(x, s.L), portal(y, s.L ), [s.M]]
        else
            x = p.r[:]
            idx = [[1], portal(x, s.L), [s.M]]
        end
        portals = vcat(idx...) |> sort |> unique
        particle_paths = []
        for (i, j) in pairwise(portals; weight = 1)
            entry = Dict()
            if s.dim == 2
                entry = Dict(:x => x[i:j], :y => y[i:j], :t => t[i:j])
            else
                entry = Dict(:x => x[i:j], :t => t[i:j])
            end
            push!(particle_paths, entry)
        end
        push!(plot_data, particle_paths)
    end
    return plot_data
end

# Below is simplified version without portals
# function extract_worldlines(s::System)::Vector{PlotData}
#     [Dict(:x => p.r[:, 1],
#           :y => p.r[:, 2],
#           :t => 1: s.M) for p in s.particles]
# end

function worldlines(s::System; linecolor = :black, dim = 2)
    d = s.dim == 2 && dim == 1 ? 1 : 2
    if d == 2
        p = plot(;fontfamily="Computer Modern", xlims = (-s.L, s.L), ylims = (-s.L, s.L), xlabel = L"x", ylabel = L"y", zlabel = L"\tau", legend = false)
    else
        p = plot(;fontfamily="Computer Modern", xlims = (-s.L, s.L), xlabel = "x", ylabel = L"\tau", legend = false)
    end
    wl_list = extract_worldlines(s)
    for wl in wl_list
        if d == 2
            for w in wl
                if all(ismissing.(w[:x])) || all(ismissing.(w[:y]))
                    continue
                end
                plot!(w[:x], w[:y], w[:t]; linecolor = linecolor, xlims=(-4.0, 4.0), ylims=(-4.0, 4.0))
            end
        else
            for w in wl
                plot!(w[:x], w[:t]; linecolor = linecolor)
            end
        end
    end
    p # The plot object can be passed to display function, i.e., display(worldlines(s))
end
function worldlines(s::System, plt::Plots.Plot{Plots.GRBackend}; linecolor = :black)
    wl_list = extract_worldlines(s)
    for wl in wl_list
        if s.dim == 2
            for w in wl
                if all(ismissing.(w[:x])) || all(ismissing.(w[:y]))
                    continue
                end
                plot!(plt, w[:x], w[:y], w[:t]; linecolor = linecolor, xlims=(-4.0, 4.0), ylims=(-4.0, 4.0))
            end
        else
            for w in wl
                plot!(plt, w[:x], w[:t]; linecolor = linecolor)
            end
        end
    end
    plt # The plot object can be passed to display function, i.e., display(worldlines(s))
end

# function worldlines(s::System; linecolor = :black, dim = :xy, xlims = (-s.L, s.L), ylims = (-s.L, s.L))
#     p = plot(legend = false)
#     wl_list = extract_worldlines(s)
#     for wl in wl_list
#         if s.dim == 2 && dim == :xy
#             plot!(p, xlims = xlims, ylims = ylims, xlabel = "x", ylabel = "y", zlabel = "τ", legend = false)
#             for w in wl
#                 if all(ismissing.(w[:x])) || all(ismissing.(w[:y]))
#                     continue
#                 end
#                 plot!(p, w[:x], w[:y], w[:t]; linecolor = linecolor)
#             end
#         elseif s.dim == 1 || dim == :x
#             plot!(p, xlims = xlims, xlabel = "x", ylabel = "τ")
#             for w in wl
#                 plot!(p, w[:x], w[:t]; linecolor = linecolor)
#             end
#         elseif dim == :y
#             plot!(p, xlims = ylims, xlabel = "y", ylabel = "τ")
#             for w in wl
#                 plot!(p, w[:y], w[:t]; linecolor = linecolor)
#             end
#         end
#     end
#     return p # The plot object can be passed to display function, i.e., display(worldlines(s))
# end

function topview(s::System)
    p = plot(;fontfamily = "Computer Modern", xlims = (-s.L, s.L), ylim = (-s.L, s.L), xlabel = L"x", ylabel = L"y", legend = false, aspect_ratio = :equal)
    for wl in extract_worldlines(s)
        for w in wl
            plot!(w[:x], w[:y]; linecolor = :black)
        end
    end
    p # Again export the plot
end

function plot_density(s::System, d::Density; readfile = false)
    if s.dim == 1
        if readfile
            df = DataFrame(CSV.File("density.csv"))
        else
            pos = [-s.L + i * d.bin + 0.5 * d.bin for i in 1:(d.nbins-1)]
            den = [d.dens[i] / (d.bin * d.ndata) for i in 1:(d.nbins-1)]
            df = DataFrame(:position => pos, :density => den)
        end

        plt = plot(df[:, 1], df[:, 2], labels = "Density", fontfamily = "Computer Modern")
        # plot!(df[:,1], s.potential, labels = "Potential")
        plot!(plt, legend = :outerright, legendfontsize = 9, foreground_color_legend = nothing, background_color_legend = nothing, ylim = (0, 5))
        plot!(plt, annotations = [
                ((1.07, 0.85), text("N = $(s.N)", :left, 9, "Computer Modern")),
                ((1.07, 0.80), text("M = $(s.M)", :left, 9, "Computer Modern")),
                ((1.07, 0.75), text("T = $(1/s.β)", :left, 9, "Computer Modern")),
                ((1.07, 0.70), text(Dates.format(now(), "dd/mm/yyyy"), :left, 9, "Computer Modern")),
                ((1.07, 0.65), text(Dates.format(now(), "HH:MM:SS"), :left, 9, "Computer Modern"))
            ], show = true)
        display(plt)
    end
    if s.dim == 2
        if readfile
            df = DataFrame(CSV.File("density.csv"))
        else
            df = density_df(s, d)
        end
        plt = heatmap(df[:, 1], df[:, 1], Matrix(df[:, Not(1)]), show = true, legend = true)
        display(plt)
    end
end

#TODO readfile for dim=2
# function plot_density(s::System, d::Density; readfile = false)
#     if s.dim == 1
#         if readfile
#             df = DataFrame(CSV.File("density.csv"))
#         else
#             pos = [-s.L + i * d.bin + 0.5 * d.bin for i in 1:(d.nbins-1)]
#             den = [d.dens[i] / (d.bin * d.ndata) for i in 1:(d.nbins-1)]
#             df = DataFrame(:position => pos, :density => den)
#         end

#         plt = plot(df[:, 1], df[:, 2], labels = "Density", fontfamily = "Computer Modern")
#         # plot!(df[:,1], s.potential, labels = "Potential")
#         plot!(plt, legend = :outerright, legendfontsize = 9, foreground_color_legend = nothing, background_color_legend = nothing, ylim = (0, 5))
#         plot!(plt, annotations = [
#                 ((1.07, 0.85), text("N = $(s.N)", :left, 9, "Computer Modern")),
#                 ((1.07, 0.80), text("M = $(s.M)", :left, 9, "Computer Modern")),
#                 ((1.07, 0.75), text("T = $(1/s.β)", :left, 9, "Computer Modern")),
#                 ((1.07, 0.70), text(Dates.format(now(), "dd/mm/yyyy"), :left, 9, "Computer Modern")),
#                 ((1.07, 0.65), text(Dates.format(now(), "HH:MM:SS"), :left, 9, "Computer Modern"))
#             ], show = true)
#         display(plt)
#     end
#     if s.dim == 2
#         if readfile
#             df = DataFrame(CSV.File("density.csv"))
#         else
#             df = density_df(s, d)
#         end
#         plt = heatmap(df[:, 1], df[:, 1], Matrix(df[:, Not(1)]), show = true, legend = true)
#         display(plt)
#     end
# end

function plot_paths(s::System; readfile = false)
    if readfile
        df = DataFrame(CSV.File("paths.csv"))
    else
        df = paths_df(s)
    end
    plt = plot(size = (800, 800))
    if s.dim == 1
        for i in 2:s.N+1
            plot!(plt, df[:, i], df[:, 1], labels = "Particle $(i-1)", fontfamily = "Computer Modern",
                legend = :outerright, foreground_color_legend = nothing, background_color_legend = nothing,
                markershape = :circle, markersize = 2, xlim = (-s.L, s.L)
            )
        end
        plot!(plt, annotations = [((1.07, 0.85), text("N = $(s.N)", :left, 9, "Computer Modern")),
                ((1.07, 0.80), text("M = $(s.M)", :left, 9, "Computer Modern")),
                ((1.07, 0.75), text("T = $(1/s.β)", :left, 9, "Computer Modern")),
                ((1.07, 0.70), text(Dates.format(now(), "dd/mm/yyyy"), :left, 9, "Computer Modern")),
                ((1.07, 0.65), text(Dates.format(now(), "HH:MM:SS"), :left, 9, "Computer Modern"))
                # ((1.07,0.60), text(join(s.permutations, " "), :left, 9, "Computer Modern"))
            ], show = true)
    elseif s.dim == 2
        for i in 2:s.N+1
            plot!(plt, df[:, 2*(i-1)], df[:, 2*(i-1)+1], df[:, 1], labels = "Particle $(i-1)",
                legend = :outerright, foreground_color_legend = nothing, background_color_legend = nothing,
                markershape = :circle, markersize = 2)
        end
    end
    display(plt)
    savefig(plt, "paths.pdf")
end
