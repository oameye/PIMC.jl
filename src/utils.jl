
unzip(a) = (getfield.(a, x) for x in fieldnames(eltype(a)))

allowmissing(x::AbstractArray{T}) where {T} = convert(AbstractArray{Union{T,Missing}}, x)
disallowmissing(x::AbstractArray{T}) where {T} = convert(AbstractArray{nonmissingtype(T)}, x)
