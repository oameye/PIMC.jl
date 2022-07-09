using ForwardDiff: derivative

function z(β, d)
    return (exp(-0.5*β)/(1-exp(-β)))^d
end;

function Z(N, β, d)
    if N == 0
        return 1.0
    else
        return (1 / N) * sum(z(k * β, d) * Z(N - k, β, d) for k in 1:N)
    end
end;

function energy(N::Int64, β::Float64, d::Int64; type::Symbol = :boson)
    if type == :boson
        return -derivative(x -> log(Z(N, x, d)), β)
    else type == :boltzmannon
        return -derivative(x -> log(z(x, d)^N), β)
    end
end;