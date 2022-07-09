using Distributed
addprocs(10)

@everywhere begin
using Pimc
using JLD2, FileIO

str(x::Float64) = replace(string(x), "."=>"_")
τ = 0.01
L = 4.0
maxval = round(√(L^2+L^2), RoundUp)
end

E = pmap(0.16:1.0:20.16, on_error=identity) do g
    store = "Store/Interactionpropagators/tau$(str(τ))/propint-L$(str(maxval))g$(str(g))tau$(str(τ)).jld2"
    println(store)
    prop_int = build_prop_int(maxval, g, τ, 1/(π^2))
    @save store prop_int
end