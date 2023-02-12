using Logging
# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debuglogger)


function info(updates::Updates)
    for (i, (n, u)) in pairs(updates)
        name = split(String(Symbol(u)), "(")[1]
        acc = acceptance(u.counter_var.queue)
        if occursin("Permutation", name)
            @info "Acceptance $(name): $(acc)"
        elseif occursin("Dummy", name)
            @info "$(name)($n):"
        elseif occursin("Bissection", name)
            @info "$(name)($n):\n\tAcceptance $(round(acc, digits=3))\n\tSlices\t   $(2^upd.level.l)"
        elseif occursin("Close", name) || occursin("Remove", name)
            @info "$(name)($n):\n\tAcceptance $(round(acc, digits=3))\n\tThreshold  $(u.var.m)"
        elseif occursin("Center", name) || occursin("Single", name)
            @info "$(name)($n):\n\tAcceptance $(round(acc, digits=3))\n\tStep\t   $(round(u.var.size, digits=3))"
        else
            @info "$(name)($n):\n\tAcceptance $(round(acc, digits=3))\n\tSlices\t   $(u.var.m)"
        end

        if i==length(updates)
            println(" ")
        end
    end
end

function info(s::System)
    # worms = s.worm_algorithm ? "Worm Algorithm" : "No Worm Algorithm"
    worms = "No Worm Algorithm"
    # ens = s.gc ? "Grand canonical" : "Canonical"
    ens = "Canonical"
    items = join([
            "Number of time slices M = $(s.M)",
            "Chemical potential    μ = $(s.μ)",
            "scattering length     a = $(s.a)",
            "Inverse temperature   β = $(s.β)",
            "Time slice size       τ = $(s.τ)",
            "$(ens) ensemble with $(worms)"], "\n")
    @info "PIMC: $(s.N) particles in a $(s.dim)D potential\n$(items)"
    println("")
end
