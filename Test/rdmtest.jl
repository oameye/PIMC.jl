using Pimc

s = System(
    _ -> 0.0;
    dim = 2,
    M = 100,
    N = 1,
    L = 10.0,
    T = 1.0,
    worms = true
)

n = 10_000_000

IR = Threshold(10, 1, s.M-1, 0.4, 0.6)
OC = Threshold(10, 1, s.M-1, 0.4, 0.6)
worms = [
    (1,  Insert(IR, adj = 1, range = 10_000)),
    (1,  Remove(IR, adj = 1, range = 10_000)),
    (1,    Open(OC, adj = 1, range = 10_000)),
    (1,   Close(OC, adj = 1, range = 10_000))
]

run!(s, n, worms)

