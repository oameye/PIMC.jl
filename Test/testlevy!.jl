using Pimc, Plots

function create_polymer(n::Int64, perm::Vector)
    pol = zeros(Int64, length(perm)) # list of particles in the polymer
    pol[1] = n
    Npol = 1
    if perm[n] != n  # go to permuations cycle
        jj = n
        while true            
            jj = perm[jj]
            if jj == n; break; end
            Npol += 1
            pol[Npol] = jj
        end
    end
    return Npol, pol
end
struct Polymer
    pol::Vector
    Npol::Int64
    n::Int64
    j₀::Int64
    function Polymer(n::Int64, perm::Vector, j₀=0)
        Npol, pol = create_polymer(n, perm)
        new(pol, Npol, n, j₀)
    end
end

function periodicbound(x::Float64, b::Float64)
    sign(x) * (mod(abs(x) + b, 2 * b) - b)
end

function perposition(b1::Vector{Float64}, b2::Vector{Float64}, xold::Vector{Float64}, size::Float64)
    x10 = xold - b1; x02 = b2 - xold
    x10 -= (2*size)*round.(x10/(2*size)); x02 -= (2*size)*round.(x02/(2*size))
    return xold - x10, xold + x02
end

function partpol(j::Int64, p::Polymer, s::System)
    return p.pol[mod1(1+floor(Int64, (j-1)/s.M), p.Npol)]
end

function partpol(j::Int64, pol, Npol, M)
    return pol[mod1(1+floor(Int64, (j-1)/M), Npol)]
end

function levy!(r′::Matrix, p::Polymer, s::System)
    m= size(r′, 2)-2
    for i in 1:m
        αᵢ = (m+1-i)/(m+2-i)
        σᵢ = √(αᵢ*s.τ)
        rᵢ, rₘ = perposition(r′[:,i], r′[:,end], s.p[partpol(p.j₀+i, p, s)].r[:,mod1(p.j₀+i, s.M)], s.size) 
        r′[:,i+1] = periodicbound.( αᵢ.*rᵢ+ (1-αᵢ).*rₘ  + randn(s.dim).*σᵢ, s.size)
    end
end

function main()
scale = 0.5
depth = 8.

function v1d(x::Vector{Float64})
    a = @. depth * sin(2π * x * scale)^2 
    return a[1]
end
s = System(v1d; dim=1, M=30, N=1, size=4., T= 1. )
yrange = [i for i in 1:s.M]
plt = plot(s.p[1].r[1, :], yrange, show=true, xlim=[-s.size ,s.size],  markershape = :circle, markersize =2)
display(plt)

n = 1; j₀ = 5
p = Polymer(n, s.permutations, j₀)
m=20
jₘ = j₀ + m

r′ = zeros(Float64, s.dim, m+1)
r′[:,1] = s.p[n].r[:, j₀]
r′[:,m+1] = s.p[partpol(jₘ, p, s)].r[:, mod1(jₘ, s.M)]
levy!(r′, p, s)
kk=[j₀+i for i in 0:m] 
plot!(plt, r′[1, :],kk, yrange, show=true, xlim=[-s.size ,s.size],  markershape = :circle, markersize =2)
display(plt)
end
main()