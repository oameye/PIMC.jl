acceptance(q::Queue{Bool})::Float64 = sum(q)/length(q)

function queue!(c::Counter, acc::Bool)::Nothing
    c.tries += 1
    enqueue!(c.queue, acc ? 1 : 0)
    if length(c.queue) > c.range
        _ = dequeue!(c.queue)
    end
    nothing
end

function apply!(s::System, fs::Vector{Measurement})::Nothing
    for f in fs
        f(s)
    end
    nothing
end
apply!(s::System, fs::ZMeasurements) = apply!(s, Measurements(fs))
function apply!(s::System, f::UpdateC)::Nothing
    acc = f(s)
    queue!(f.counter, acc)

    if f.counter_var.tries % f.counter_var.adj == 0
        adjust!(f.var, acceptance(f.counter_var.queue))
    end
    nothing
end

function run!(s::System, n::Int64, updates::Updates; Zmeasurements::ZMeasurements = ZMeasurement[])::Nothing
    # number of tries to sample beads outside hardsphere; smaller than in thermalization phase
    therm = isempty(Zmeasurements) ? true : false
    s.ctr = therm ? 10_000 : 1_000
    (every, fs) = unzip(updates)
    weights = Weights(1 ./ every)
    for _ = 1:n
        f = sample(fs, weights)
        apply!(s, f)

        measurement_Z_sector(s, Zmeasurements)
    end
    nothing
end
