@testset "Measurements" begin

@testset "Density" begin
	@testset "Number of particles with potential" begin
        scale = 0.5
        depth = 8.

        v1d(x::Vector{Float64}) = depth * sin(2π * x[1] * scale)^2

		s = Pimc.System(v1d)

		updates = [
			(2, SingleCenterOfMass(s, 3.0)),
			(1, ReshapeLinear(s, 20)),
			(1, ReshapeSwapLinear(s, 20))
		]

        d = Density(s)
        measurements = ZMeasurement[d]

		run!(s, 10000, updates, Zmeasurements=measurements)

		# sanity check
		numb = (sum(d.dens) / d.ndata)
		if abs(numb - s.N) > 1e-2
			@error "Error save_density: integrated density = $numb, should be $(s.N) if x-range is ok"
		end

		@test isapprox(numb, s.N, atol=1e-2)
	end
end

end # @testset "Measurements"
