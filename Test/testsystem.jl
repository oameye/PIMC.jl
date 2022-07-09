# if 2*s.L/s.nbins < s.a
#     @info "binsize ($(2*s.L/s.nbins)) is smaller than 2d scattering length ($(s.a))"
#     return
# end

function normalized_intensity(coord::Vector{Float64}, angles::Vector{Float64}, scale::Float64)::Float64
    s, c = 0.0, 0.0

    # Compute plane wave contributions for all laser beam orientations
    for angle in angles
        r = coord[1] * sin(angle) + coord[2] * cos(angle)
        s += sin(2π * r * scale)
        c += cos(2π * r * scale)
    end

    # Normalize max intensity to be 1
    s /= length(angles)
    c /= length(angles)

    return (s * s + c * c)
end

@testset "System" begin

@testset "Optical Lattice Potential" begin
	@testset "Calculation" begin
		# Initial conditions with known solution
		cyclic_3 = [2π*k/3 for k in 0:2]
		scale = 1.
		depth = 1.

		# For testing numerics use \approx
		@test normalized_intensity([0.0, 0.0], cyclic_3, scale) ≈ 1.0

		# Test 3 peaks of the unit cell
		@test normalized_intensity([2/√3, 0.0], cyclic_3, scale) ≈ 1.0
		@test normalized_intensity([0.0, 2/3], cyclic_3, scale) ≈ 1.0
	end
end

@testset "System setup" begin
	@testset "Initialization1d" begin
        scale = 0.5
        depth = 8.

        v1d(x::Vector{Float64}) = depth * sin(2π * x[1] * scale)^2

		s = Pimc.System(v1d)

		# Dummy test
		@test s.N == 2
	end

    @testset "Initialization2d" begin
        function normalized_intensity(coord::Vector{Float64}, angles::Vector{Float64}, scale::Float64)::Float64
            s, c = 0.0, 0.0
        
            # Compute plane wave contributions for all laser beam orientations
            for angle in angles
                r = coord[1] * sin(angle) + coord[2] * cos(angle)
                s += sin(2π * r * scale)
                c += cos(2π * r * scale)
            end
        
            # Normalize max intensity to be 1
            s /= length(angles)
            c /= length(angles)
        
            return (s * s + c * c)
        end
        
        function generate_V(angles::Vector{Float64}, scale::Float64, depth::Float64)::Function
            V(coord::Vector{Float64}) = depth*normalized_intensity(coord, angles, scale)
        end

        cyclic_4 = [2π*k/4 for k in 0:3]
        scale = 0.5
        depth = 8.
        potential = generate_V(cyclic_4, scale, depth)

        s = System(potential; dim=2, M=100, N=5, L=4., T=1)

		# Dummy test
		@test s.N == 5
	end
end

@testset "Periodic Bounds" begin
    s = Pimc.System(x -> 0)
    updates = [
        (2, CenterOfMass(s, 3.0)),
        (1, ReshapeLinear(s, 20)),
        (1, ReshapeSwapLinear(s, 20))
    ]
    
    run!(s, 10_000, updates)

    check = 0
    for n in 1:s.N
        for m in 1:s.M
            if all(s.world[n].r[m, :] .> s.L) ||  all(s.world[n].r[m, :] .< -s.L)
                check+=1
            end
        end
    end
    @test check == 0
end

end # @testset "system"