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

@testset "Optical Lattice Potential" begin
	@testset "Calculation" begin
		# Initial conditions with known solution
		cyclic_3 = [2π*k/3 for k in 0:2]
		scale = 1.
		depth = 1.

		# For testing numerics use \approx (≈)
		@test normalized_intensity([0.0, 0.0], cyclic_3, scale) ≈ 1.0

		# Test 3 peaks of the unit cell
		@test normalized_intensity([2/√3, 0.0], cyclic_3, scale) ≈ 1.0
		@test normalized_intensity([0.0, 2/3], cyclic_3, scale) ≈ 1.0
	end
end
