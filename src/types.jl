abstract type Worldline end
abstract type Measurement end
abstract type ZMeasurement <: Measurement end
abstract type Update end
abstract type UpdateC <: Update end

# Perform updates/measurements to the system.
# Each function comes with an integer that denotes the frequency of the Update.
Updates = Vector{Tuple{Int64,A}} where {A <: Update}
ZMeasurements = Vector{ZMeasurement}
Measurements = Vector{Measurement}

VectorMissing{T} = Vector{Union{Missing,T}}
