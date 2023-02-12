using Pkg
current_path = @__DIR__
Pkg.activate(current_path * "/../.");

using Pimc, Test, Logging

files = [
    "testpotential.jl",
    "testsystem.jl",
    "testmeasurements.jl",
    "testupdates.jl"
    ]

for file in files
    include(file)
    # printstyled(file * ":    OK\n"; color = :green)
end

printstyled("\nALL TESTS PASSED!\n"; color = :green)
