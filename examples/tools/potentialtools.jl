function normalized_intensity(coord::Vector{Float64}, angles::Vector{Float64}, scale::Float64; helical::Bool = false)::Float64
    s, c = 0.0, 0.0

    # Compute plane wave contributions for all laser beam orientations
    for angle in angles
        r = coord[1] * sin(angle) + coord[2] * cos(angle)
        s += helical ? sin(2π * r * scale + angle) : sin(2π * r * scale)
        c += helical ? cos(2π * r * scale + angle) : cos(2π * r * scale)
    end

    # Normalize max intensity to be 1
    s /= length(angles)
    c /= length(angles)

    return (s * s + c * c)
end

function generate_V(angles::Vector{Float64}, scale::Float64, depth::Float64; sym::Symbol=:sim, helical::Bool = false)::Function
    if sym ==:plot
        return V(coord1::Float64, coord2::Float64) = depth * normalized_intensity([coord1, coord2], angles, scale, helical=helical)
    else
        return V(coord::Vector{Float64}) = depth * normalized_intensity(coord, angles, scale, helical)
    end
end
function generate_V(scale::Float64, depth::Float64, sym::Symbol; attractive = true)::Function
    if sym==:l65
        ang = [0.5191461142465229, 1.695151321341658, -0.12435499454676144, -1.4464413322481353, 2.0899424410414196, 2.62244653934327, 0.12435499454676144, -1.0516502125483738, -1.695151321341658, 1.446441332248135, -2.0899424410414196, -3.017237659043032, -0.5191461142465229, 3.017237659043032, -2.62244653934327, 1.0516502125483738]
    elseif sym == :l25
        ang = [2.214297435588181,0.9272952180016122,-0.6435011087932844,0.6435011087932844,-2.498091544796509,3.141592653589793, 2.498091544796509,0,1.5707963267948966,-2.2142974355881813,-1.5707963267948968,-0.9272952180016123]
    elseif sym == :cubic
        ang = [2π * k / 4 for k in 0:3]
    elseif sym == :harmonic
        return (r) -> 0.5*(r[1]^2+r[2]^2)
    else
        error()
    end
    sgn = attractive ? -1 : 1
    return V(coord::Vector{Float64}) = sgn*depth * normalized_intensity(coord, ang, scale)
end

# function meshgrid(x::LinRange{Float64, Int64}, y::LinRange{Float64, Int64})::Tuple{Matrix{Float64}, Matrix{Float64}}
#     X = [x for _ in y, x in x]
#     Y = [y for y in y, _ in x]
#     X, Y
# end

function meshgrid(x, y)
    X = [x for _ in y, x in x]
    Y = [y for y in y, _ in x]
    X, Y
end


l25 = [2.214297435588181,
0.9272952180016122,
-0.6435011087932844,
0.6435011087932844,
-2.498091544796509,
3.141592653589793,
2.498091544796509,
0.0,
1.5707963267948966,
-2.2142974355881813,
-1.5707963267948968,
-0.9272952180016123]
l65 = [0.5191461142465229,
    1.695151321341658,
    -0.12435499454676144,
    -1.4464413322481353,
    2.0899424410414196,
    2.62244653934327,
    0.12435499454676144,
    -1.0516502125483738,
    -1.695151321341658,
    1.446441332248135,
    -2.0899424410414196,
    -3.017237659043032,
    -0.5191461142465229,
    3.017237659043032,
    -2.62244653934327,
    1.0516502125483738
]