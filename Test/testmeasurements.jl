@testset "Measurements" begin

@testset "Density" begin
	@testset "Number of particles with potential" begin
        scale = 0.5
        depth = 8.

        v1d(x::Vector{Float64}) = depth * sin(2π * x[1] * scale)^2

		s = Pimc.System(v1d)

		updates = [
			(2, CenterOfMass(s, 3.0)),
			(1, ReshapeLinear(s, 20)),
			(1, ReshapeSwapLinear(s, 20))
		]

		d = Density(s)
		measurements = [
    		(10, d)
			]

		run!(s, 10000, updates, measurements)

		# sanity check
		numb = (sum(d.dens) / d.ndata)
		if abs(numb - s.N) > 1e-4
			@error "Error save_density: integrated density = $numb, should be $(s.N) if x-range is ok"
		end

		@test numb ≈ s.N atol=1e-1
	end
end

end # @testset "Measurements" 