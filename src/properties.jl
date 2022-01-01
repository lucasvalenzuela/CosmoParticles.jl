"""
    colnorm(a::AbstractMatrix[, origin::AbstractVector])

Returns a new `Vector` with the columnwise norms of `a` around `origin`.

This method does not check for overflow of the squared norm.
"""
function colnorm(a::AbstractMatrix{T}) where {T}
    dst = Vector{T}(undef, size(a, 2))
    return _colnorm_unsafe!(dst, a)
end

function colnorm(a::AbstractMatrix{T}, origin::AbstractVector{T0}) where {T,T0}
    @assert size(a, 1) == length(origin)
    T === T0 && (origin = convert(Vector{T}, origin))
    dst = Vector{T}(undef, size(a, 2))
    return _colnorm_unsafe!(dst, a, origin)
end

"""
    colnorm2(a::AbstractMatrix[, origin:AbstractVector])

Returns a new `Vector` with the columnwise squared norms of `a` around `origin`.
"""
function colnorm2(a::AbstractMatrix{T}) where {T}
    dst = Vector{T}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a)
end

function colnorm2(a::AbstractMatrix{Q}) where {T,Q<:Quantity{T}}
    u = unit(Q)^2
    Q² = Quantity{T,dimension(u),typeof(u)}
    dst = Vector{Q²}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a)
end

function colnorm2(a::AbstractMatrix{Union{Missing,Q}}) where {T,Q<:Quantity{T}}
    u = unit(Q)^2
    Q² = Quantity{T,dimension(u),typeof(u)}
    dst = Vector{Union{Missing,Q²}}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a)
end

function colnorm2(a::AbstractMatrix{T}, origin::AbstractVector{T0}) where {T,T0}
    @assert size(a, 1) == length(origin)
    T === T0 && (origin = convert(Vector{T}, origin))
    dst = Vector{T}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, origin)
end

function colnorm2(a::AbstractMatrix{Q}, origin::AbstractVector{Q0}) where {T,Q<:Quantity{T},Q0}
    @assert size(a, 1) == length(origin)
    Q === Q0 && (origin = convert(Vector{Q}, origin))
    u = unit(Q)^2
    Q² = Quantity{T,dimension(u),typeof(u)}
    dst = Vector{Q²}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, origin)
end

function colnorm2(a::AbstractMatrix{Union{Missing,Q}}, origin::AbstractVector{Q0}) where {T,Q<:Quantity{T},Q0}
    @assert size(a, 1) == length(origin)
    Q === Q0 && (origin = convert(Vector{Q}, origin))
    u = unit(Q)^2
    Q² = Quantity{T,dimension(u),typeof(u)}
    dst = Vector{Union{Missing,Q²}}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, origin)
end


function _colnorm_unsafe! end
function _colnorm2_unsafe! end

for (name, zeroT2) in zip(["_colnorm_unsafe", "_colnorm2_unsafe"], [:(zero(T)^2), :(zero(T))])
    quote
        function $(Symbol(name, "!"))(dst::AbstractVector{T}, a::AbstractMatrix) where {T}
            dims = size(a, 1)
            if dims == 2
                return $(Symbol(name, "_2D!"))(dst, a)
            elseif dims == 3
                return $(Symbol(name, "_3D!"))(dst, a)
            else
                return $(Symbol(name, "_ND!"))(dst, a; zeroT2=$zeroT2)
            end
        end

        function $(Symbol(name, "!"))(
            dst::AbstractVector{T},
            a::AbstractMatrix,
            origin::AbstractVector,
        ) where {T}
            dims = size(a, 1)
            if dims == 2
                return $(Symbol(name, "_2D!"))(dst, a, origin)
            elseif dims == 3
                return $(Symbol(name, "_3D!"))(dst, a, origin)
            else
                return $(Symbol(name, "_ND!"))(dst, a, origin; zeroT2=$zeroT2)
            end
        end
    end |> eval
end

function _colnorm2_unsafe_ND!(dst::AbstractVector{T}, a::AbstractMatrix; zeroT2=zero(T)^2) where {T}
    @inbounds @simd for i in axes(a, 2)
        dsti = zeroT2
        for j in axes(a, 1)
            dsti += a[j, i]^2
        end
        dst[i] = dsti
    end
    return dst
end

function _colnorm2_unsafe_2D!(dst::AbstractVector, a::AbstractMatrix)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = a[1, i]^2 + a[2, i]^2
    end
    return dst
end

function _colnorm2_unsafe_3D!(dst::AbstractVector, a::AbstractMatrix)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = a[1, i]^2 + a[2, i]^2 + a[3, i]^2
    end
    return dst
end

function _colnorm2_unsafe_ND!(
    dst::AbstractVector{T},
    a::AbstractMatrix,
    origin::AbstractVector;
    zeroT2=zero(T)^2,
) where {T}
    @inbounds @simd for i in axes(a, 2)
        dsti = zeroT2
        for j in axes(a, 1)
            dsti += (a[j, i] - origin[j])^2
        end
        dst[i] = dsti
    end
    return dst
end

function _colnorm2_unsafe_2D!(dst::AbstractVector, a::AbstractMatrix, origin::AbstractVector)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = (a[1, i] - origin[1])^2 + (a[2, i] - origin[2])^2
    end
    return dst
end

function _colnorm2_unsafe_3D!(dst::AbstractVector, a::AbstractMatrix, origin::AbstractVector)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = (a[1, i] - origin[1])^2 + (a[2, i] - origin[2])^2 + (a[3, i] - origin[3])^2
    end
    return dst
end

function _colnorm_unsafe_ND!(dst::AbstractVector{T}, a::AbstractMatrix; zeroT2=zero(T)^2) where {T}
    @inbounds @simd for i in axes(a, 2)
        dsti = zeroT2
        for j in axes(a, 1)
            dsti += a[j, i]^2
        end
        dst[i] = sqrt(dsti)
    end
    return dst
end

function _colnorm_unsafe_2D!(dst::AbstractVector, a::AbstractMatrix)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = sqrt(a[1, i]^2 + a[2, i]^2)
    end
    return dst
end

function _colnorm_unsafe_3D!(dst::AbstractVector, a::AbstractMatrix)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = sqrt(a[1, i]^2 + a[2, i]^2 + a[3, i]^2)
    end
    return dst
end

function _colnorm_unsafe_ND!(
    dst::AbstractVector{T},
    a::AbstractMatrix,
    origin::AbstractVector;
    zeroT2=zero(T)^2,
) where {T}
    @inbounds @simd for i in axes(a, 2)
        dsti = zeroT2
        for j in axes(a, 1)
            dsti += (a[j, i] - origin[j])^2
        end
        dst[i] = sqrt(dsti)
    end
    return dst
end

function _colnorm_unsafe_2D!(dst::AbstractVector, a::AbstractMatrix, origin::AbstractVector)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = sqrt((a[1, i] - origin[1])^2 + (a[2, i] - origin[2])^2)
    end
    return dst
end

function _colnorm_unsafe_3D!(dst::AbstractVector, a::AbstractMatrix, origin::AbstractVector)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = sqrt((a[1, i] - origin[1])^2 + (a[2, i] - origin[2])^2 + (a[3, i] - origin[3])^2)
    end
    return dst
end
