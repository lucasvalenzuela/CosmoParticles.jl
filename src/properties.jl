"""
    meanprop(p::AbstractParticles, prop::Symbol; massprop=:mass, massweighted::Bool=true)

Returns the mean of the specified property over all particles.

The mean can be weighted by the mass.
"""
function meanprop(p::AbstractParticles, prop::Symbol; massprop=:mass, massweighted=true)
    if massweighted
        return _meanprop(p[prop], p[massprop])
    else
        return _meanprop(p[prop])
    end
end

_meanprop(a::Number) = a
_meanprop(a::Number, ::Any) = a
_meanprop(a::Fill) = a[1]
_meanprop(a::Fill, ::AbstractVector) = a[1]
_meanprop(a::Fill, ::Fill) = a[1]
_meanprop(a::AbstractVector) = mean(a)
_meanprop(a::AbstractMatrix) = mean(a; dims=2)
_meanprop(a::AbstractVector, m::AbstractVector) = mean(a, weights(ustrip_lazy(m)))
_meanprop(a::AbstractMatrix, m::AbstractVector) = mean(a, weights(ustrip_lazy(m)); dims=2)
_meanprop(a::AbstractArray, ::Number) = _meanprop(a)
_meanprop(a::AbstractVector, ::Fill) = _meanprop(a)
_meanprop(a::AbstractMatrix, ::Fill) = _meanprop(a)


"""
    sumprop(p::AbstractParticles, prop::Symbol)

Returns the sum of the specified property over all particles.
"""
sumprop(p::AbstractParticles, prop::Symbol) = _sumprop(p, p[prop])

_sumprop(p::AbstractParticles, a::Number) = particle_number(p) * a
_sumprop(::AbstractParticles, a::AbstractVecOrMat) = _sumprop(a)
_sumprop(a::AbstractVector) = sum(a)
_sumprop(a::AbstractMatrix) = sum(a; dims=2)


"""
    meanpos(p::AbstractParticles, prop::Symbol=:pos; massprop=:mass, massweighted=true)

Returns the mean of the particles' positions.

The mean can be weighted by the mass.
"""
meanpos(p::AbstractParticles, prop::Symbol=:pos; massprop=:mass, massweighted=true) =
    meanprop(p, prop; massprop, massweighted)

"""
    meanvel(p::AbstractParticles, prop::Symbol=:pos; massprop=:mass, massweighted=true)

Returns the mean of the particles' velocities.

The mean can be weighted by the mass.
"""
meanvel(p::AbstractParticles, prop::Symbol=:vel; massprop=:mass, massweighted=true) =
    meanprop(p, prop; massprop, massweighted)


"""
    angmom(p::AbstractParticles)

Returns the individual angular momenta of the particles.

If `angmomprop` is a `Symbol`, the already computed angular momenta are returned directly if they exist.

# Keyword Arguments
- `origin=Nothing`: origin of the coordinate system for the computation of the angular momentum
- `velorigin=Nothing`: relative velocity zero-point, oftentimes the velocity of the center of mass
- `posprop=:pos`: property name of the position
- `velprop=:vel`: property name of the velocity
- `massprop=:mass`: property name of the mass
- `angmomprop=nothing`: property name of the angular momentum - the method checks whether it already exists or not. Warning: The method assumes the already saved property to have used the same `origin` and `velorigin`.
"""
function angmom(
    p::AbstractParticles;
    origin=nothing,
    velorigin=nothing,
    posprop=:pos,
    velprop=:vel,
    massprop=:mass,
    angmomprop=nothing,
)
    !isnothing(angmomprop) && haskey(p, angmomprop) && return p[angmomprop]
    return angmom(p[posprop], p[velprop], p[massprop]; origin, velorigin)
end

"""
    angmomtot(p::AbstractParticles)

Returns the total angular momentum of the particles.

If `angmomprop` is a `Symbol`, the total angular momentum is computed from the already existing angular momenta.

# Keyword Arguments
- `origin=Nothing`: origin of the coordinate system for the computation of the angular momentum
- `velorigin=Nothing`: relative velocity zero-point, oftentimes the velocity of the center of mass
- `posprop=:pos`: property name of the position
- `velprop=:vel`: property name of the velocity
- `massprop=:mass`: property name of the mass
- `angmomprop=nothing`: property name of the angular momentum - the method checks whether it already exists or not. Warning: The method assumes the already saved property to have used the same `origin` and `velorigin`.
"""
function angmomtot(
    p::AbstractParticles;
    origin=nothing,
    velorigin=nothing,
    posprop=:pos,
    velprop=:vel,
    massprop=:mass,
    angmomprop=nothing,
)
    !isnothing(angmomprop) && haskey(p, angmomprop) && return sumprop(p, angmomprop)
    return angmomtot(p[posprop], p[velprop], p[massprop]; origin, velorigin)
end

"""
    angmomtot_stable(p::AbstractParticles)

Returns the total angular momentum of the particles, like [`angmomtot`](@ref).

Uses the stable summation algorithm `Base.sum` for the summation, but is almost three times slower than
[`angmomtot`](@ref).

See [`angmomtot`](@ref) for the available keyword arguments.
"""
function angmomtot_stable(
    p::AbstractParticles;
    origin=nothing,
    velorigin=nothing,
    posprop=:pos,
    velprop=:vel,
    massprop=:mass,
    angmomprop=nothing,
)
    !isnothing(angmomprop) && haskey(p, angmomprop) && return sumprop(p, angmomprop)
    return angmomtot_stable(p[posprop], p[velprop], p[massprop]; origin, velorigin)
end



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



"""
    coldot(a::AbstractMatrix, b::AbstractMatrix)

Returns the dot products of the matrix columns as a new matrix.

Both input matrices need to have the same dimensions ``d × N``, representing ``N`` ``d``-dimensional vectors.
"""
function coldot(a::AbstractMatrix, b::AbstractMatrix)
    @assert size(a) == size(b)
    return _coldot_unsafe(a, b)
end

function _coldot_unsafe(a::AbstractMatrix, b::AbstractMatrix)
    dst = similar(a, Any, size(a, 2))
    return _coldot_unsafe!(dst, a, b)
end

function _coldot_unsafe(a::AbstractMatrix{T1}, b::AbstractMatrix{T2}) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2))
    dst = similar(a, T, size(a, 2))
    return _coldot_unsafe!(dst, a, b)
end

function _coldot_unsafe(
    a::AbstractMatrix{Union{T1,Missing}},
    b::AbstractMatrix{T2},
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2))
    dst = similar(a, Union{T,Missing}, size(a, 2))
    return _coldot_unsafe!(dst, a, b)
end

function _coldot_unsafe(
    a::AbstractMatrix{T1},
    b::AbstractMatrix{Union{T2,Missing}},
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2))
    dst = similar(a, Union{T,Missing}, size(a, 2))
    return _coldot_unsafe!(dst, a, b)
end

function _coldot_unsafe(
    a::AbstractMatrix{Union{T1,Missing}},
    b::AbstractMatrix{Union{T2,Missing}},
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2))
    dst = similar(a, Union{T,Missing}, size(a, 2))
    return _coldot_unsafe!(dst, a, b)
end

function _coldot_unsafe!(dst::AbstractVector, a::AbstractMatrix, b::AbstractMatrix)
    dims = size(a, 1)
    if dims == 2
        return _coldot_unsafe_2D!(dst, a, b)
    elseif dims == 3
        return _coldot_unsafe_3D!(dst, a, b)
    else
        return _coldot_unsafe_ND!(dst, a, b)
    end
end

function _coldot_unsafe_2D!(dst::AbstractVector, a::AbstractMatrix, b::AbstractMatrix)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = a[1, i] * b[1, i] + a[2, i] * b[2, i]
    end
    return dst
end

function _coldot_unsafe_3D!(dst::AbstractVector, a::AbstractMatrix, b::AbstractMatrix)
    @inbounds @simd for i in axes(a, 2)
        dst[i] = a[1, i] * b[1, i] + a[2, i] * b[2, i] + a[3, i] * b[3, i]
    end
    return dst
end

function _coldot_unsafe_ND!(dst::AbstractVector{T}, a::AbstractMatrix, b::AbstractMatrix) where {T}
    @inbounds for i in axes(a, 2)
        dsti = a[1, i] * b[1, i]
        for j in Iterators.drop(axes(a, 1), 1)
            dsti += a[j, i] * b[j, i]
        end
        dst[i] = dsti
    end
    return dst
end



"""
    colcross(a::AbstractMatrix, b::AbstractMatrix)

Returns the cross products of the matrix columns as a new matrix.

Both input matrices need to have the same dimensions ``3 × N``.
"""
function colcross(a::AbstractMatrix, b::AbstractMatrix)
    @assert size(a) == size(b)
    @assert size(a, 1) == 3
    return _colcross_unsafe(a, b)
end

function _colcross_unsafe(a::AbstractMatrix, b::AbstractMatrix)
    dst = similar(a, Any)
    return _colcross_unsafe!(dst, a, b)
end

function _colcross_unsafe(a::AbstractMatrix{T1}, b::AbstractMatrix{T2}) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2))
    dst = similar(a, T)
    return _colcross_unsafe!(dst, a, b)
end

function _colcross_unsafe(
    a::AbstractMatrix{Union{T1,Missing}},
    b::AbstractMatrix{T2},
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2))
    dst = similar(a, Union{T,Missing})
    return _colcross_unsafe!(dst, a, b)
end

function _colcross_unsafe(
    a::AbstractMatrix{T1},
    b::AbstractMatrix{Union{T2,Missing}},
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2))
    dst = similar(a, Union{T,Missing})
    return _colcross_unsafe!(dst, a, b)
end

function _colcross_unsafe(
    a::AbstractMatrix{Union{T1,Missing}},
    b::AbstractMatrix{Union{T2,Missing}},
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2))
    dst = similar(a, Union{T,Missing})
    return _colcross_unsafe!(dst, a, b)
end

function _colcross_unsafe!(dst::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix)
    @inbounds @simd for i in axes(dst, 2)
        dst[1, i] = a[2, i] * b[3, i] - a[3, i] * b[2, i]
        dst[2, i] = a[3, i] * b[1, i] - a[1, i] * b[3, i]
        dst[3, i] = a[1, i] * b[2, i] - a[2, i] * b[1, i]
    end
    return dst
end



"""
    angmom(pos::AbstractMatrix, vel::AbstractMatrix, mass; origin=nothing, velorigin=nothing)

Returns the angular momentum based on the positions, velocities, and masses.

Both input matrices need to have the same dimensions ``3 × N``.
The mass can be a vector or a scalar value (useful if all particles have the same mass).

The coordinate system's origin can be passed as a vector of length 3, as well as the velocity origin, which
determines the relative velocity at which the angular momentum should be determined. Typically, this will be the
velocity of the center of mass.
"""
function angmom(
    pos::AbstractMatrix,
    vel::AbstractMatrix,
    mass::AbstractVector;
    origin=nothing,
    velorigin=nothing,
)
    @assert size(pos, 2) == size(vel, 2) == length(mass)
    @assert size(pos, 1) == size(vel, 1) == 3
    return _angmom_unsafe(pos, vel, mass, origin, velorigin)
end

function angmom(pos::AbstractMatrix, vel::AbstractMatrix, mass::Number; origin=nothing, velorigin=nothing)
    @assert size(pos, 2) == size(vel, 2)
    @assert size(pos, 1) == size(vel, 1) == 3
    return _angmom_unsafe(pos, vel, mass, origin, velorigin)
end

function _angmom_unsafe(pos::AbstractMatrix, vel::AbstractMatrix, mass, origin, velorigin)
    dst = similar(pos, Any)
    return _angmom_unsafe!(dst, pos, vel, mass, origin, velorigin)
end


_help_convert_origin(::Type{T}, origin::AbstractVector) where {T} = convert(Vector{T}, origin)
_help_convert_origin(::Type{T}, origin::Nothing) where {T} = nothing


function _angmom_unsafe(
    pos::AbstractMatrix{T1},
    vel::AbstractMatrix{T2},
    mass,
    origin,
    velorigin,
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2) * _zero_of(mass))
    dst = similar(pos, T)
    return _angmom_unsafe!(
        dst,
        pos,
        vel,
        mass,
        _help_convert_origin(T1, origin),
        _help_convert_origin(T2, velorigin),
    )
end

function _angmom_unsafe(
    pos::AbstractMatrix{Union{T1,Missing}},
    vel::AbstractMatrix{T2},
    mass,
    origin,
    velorigin,
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2) * _zero_of(mass))
    dst = similar(pos, Union{T,Missing})
    return _angmom_unsafe!(dst, pos, vel, mass, origin, velorigin)
end

function _angmom_unsafe(
    pos::AbstractMatrix{T1},
    vel::AbstractMatrix{Union{T2,Missing}},
    mass,
    origin,
    velorigin,
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2) * _zero_of(mass))
    dst = similar(pos, Union{T,Missing})
    return _angmom_unsafe!(dst, pos, vel, mass, origin, velorigin)
end

function _angmom_unsafe(
    pos::AbstractMatrix{Union{T1,Missing}},
    vel::AbstractMatrix{Union{T2,Missing}},
    mass,
    origin,
    velorigin,
) where {T1<:Number,T2<:Number}
    T = typeof(zero(T1) * zero(T2) * _zero_of(mass))
    dst = similar(pos, Union{T,Missing})
    return _angmom_unsafe!(dst, pos, vel, mass, origin, velorigin)
end

for (T1, T2) in zip(
    [:Nothing, :AbstractVector, :Nothing, :AbstractVector],
    [:Nothing, :Nothing, :AbstractVector, :AbstractVector],
)
    for T3 in [:Number, :AbstractVector]
        if T1 === :Nothing
            p1 = :(p[1, i])
            p2 = :(p[2, i])
            p3 = :(p[3, i])
        else
            p1 = :((p[1, i] - origin[1]))
            p2 = :((p[2, i] - origin[2]))
            p3 = :((p[3, i] - origin[3]))
        end
        if T2 === :Nothing
            v1 = :(v[1, i])
            v2 = :(v[2, i])
            v3 = :(v[3, i])
        else
            v1 = :((v[1, i] - velorigin[1]))
            v2 = :((v[2, i] - velorigin[2]))
            v3 = :((v[3, i] - velorigin[3]))
        end
        if T3 === :Number
            mi = :m
        else
            mi = :(m[i])
        end

        quote
            function _angmom_unsafe!(
                dst::AbstractMatrix,
                p::AbstractMatrix,
                v::AbstractMatrix,
                m::$T3,
                origin::$T1,
                velorigin::$T2,
            )
                @inbounds @simd for i in axes(dst, 2)
                    dst[1, i] = $mi * ($p2 * $v3 - $p3 * $v2)
                    dst[2, i] = $mi * ($p3 * $v1 - $p1 * $v3)
                    dst[3, i] = $mi * ($p1 * $v2 - $p2 * $v1)
                end
                return dst
            end
        end |> eval
    end
end


_zero_of(::T) where {T<:Number} = zero(T)
_zero_of(::AbstractArray{T}) where {T<:Number} = zero(T)



"""
    angmomtot(pos::AbstractMatrix, vel::AbstractMatrix, mass; origin=nothing, velorigin=nothing)

Returns the total angular momentum based on the positions, velocities, and masses.

Both input matrices need to have the same dimensions ``3 × N``.
The mass can be a vector or a scalar value (useful if all particles have the same mass).

The coordinate system's origin can be passed as a vector of length 3, as well as the velocity origin, which
determines the relative velocity at which the angular momentum should be determined. Typically, this will be the
velocity of the center of mass.

This method never allocates the full angular momentum matrix for all particles like [`angmom`](@ref) does.
Use this when the individual particle angular momenta are not needed.
"""
function angmomtot(
    pos::AbstractMatrix,
    vel::AbstractMatrix,
    mass::AbstractVector;
    origin=nothing,
    velorigin=nothing,
)
    @assert size(pos, 2) == size(vel, 2) == length(mass)
    @assert size(pos, 1) == size(vel, 1) == 3
    return _angmomtot_unsafe(pos, vel, mass, origin, velorigin)
end

function angmomtot(pos::AbstractMatrix, vel::AbstractMatrix, mass::Number; origin=nothing, velorigin=nothing)
    @assert size(pos, 2) == size(vel, 2)
    @assert size(pos, 1) == size(vel, 1) == 3
    return _angmomtot_unsafe(pos, vel, mass, origin, velorigin)
end

function _angmomtot_unsafe(pos::AbstractMatrix, vel::AbstractMatrix, mass, origin, velorigin)
    T, zeroT = _type_and_zero_of_summed_angmomtot(pos, vel, mass)
    dst = Vector{T}(undef, 3)
    _angmomtot_unsafe!(dst, pos, vel, mass, origin, velorigin; zeroT)
end

function _type_and_zero_of_summed_angmomtot(
    p::AbstractMatrix{T1},
    v::AbstractMatrix{T2},
    m::T3,
) where {T1,T2,T3<:Number}
    z = zero(T1) * zero(T2)
    return typeof(z * zero(T3)), z
end

function _type_and_zero_of_summed_angmomtot(
    p::AbstractMatrix{T1},
    v::AbstractMatrix{T2},
    m::AbstractVector{T3},
) where {T1,T2,T3}
    z = zero(T1) * zero(T2) * zero(T3)
    return typeof(z), z
end


"""
    angmomtot_stable(pos::AbstractMatrix, vel::AbstractMatrix, mass; origin=nothing, velorigin=nothing)

Returns the total angular momentum based on the positions, velocities, and masses, like [`angmomtot`](@ref).

Uses the stable summation algorithm `Base.sum` for the summation, but is almost three times slower than
[`angmomtot`](@ref).

See [`angmomtot`](@ref) for the available keyword arguments.
"""
function angmomtot_stable end

function angmomtot_stable(
    pos::AbstractMatrix,
    vel::AbstractMatrix,
    mass::AbstractVector;
    origin=nothing,
    velorigin=nothing,
)
    @assert size(pos, 2) == size(vel, 2) == length(mass)
    @assert size(pos, 1) == size(vel, 1) == 3
    return _angmomtot_stable_unsafe(pos, vel, mass, origin, velorigin)
end

function angmomtot_stable(
    pos::AbstractMatrix,
    vel::AbstractMatrix,
    mass::Number;
    origin=nothing,
    velorigin=nothing,
)
    @assert size(pos, 2) == size(vel, 2)
    @assert size(pos, 1) == size(vel, 1) == 3
    return _angmomtot_stable_unsafe(pos, vel, mass, origin, velorigin)
end



for (T1, T2) in zip(
    [:Nothing, :AbstractVector, :Nothing, :AbstractVector],
    [:Nothing, :Nothing, :AbstractVector, :AbstractVector],
)
    if T1 === :Nothing
        p1 = :(p[1, i])
        p2 = :(p[2, i])
        p3 = :(p[3, i])
    else
        p1 = :((p[1, i] - origin[1]))
        p2 = :((p[2, i] - origin[2]))
        p3 = :((p[3, i] - origin[3]))
    end
    if T2 === :Nothing
        v1 = :(v[1, i])
        v2 = :(v[2, i])
        v3 = :(v[3, i])
    else
        v1 = :((v[1, i] - velorigin[1]))
        v2 = :((v[2, i] - velorigin[2]))
        v3 = :((v[3, i] - velorigin[3]))
    end

    quote
        function _angmomtot_unsafe!(
            dst::AbstractVector,
            p::AbstractMatrix,
            v::AbstractMatrix,
            m::Number,
            origin::$T1,
            velorigin::$T2;
            zeroT,
        )
            d1 = d2 = d3 = zeroT
            @inbounds @simd for i in axes(p, 2)
                d1 += $p2 * $v3 - $p3 * $v2
                d2 += $p3 * $v1 - $p1 * $v3
                d3 += $p1 * $v2 - $p2 * $v1
            end
            dst .= m * d1, m * d2, m * d3
            return dst
        end

        function _angmomtot_unsafe!(
            dst::AbstractVector,
            p::AbstractMatrix,
            v::AbstractMatrix,
            m::AbstractVector,
            origin::$T1,
            velorigin::$T2;
            zeroT,
        )
            dst .= zeroT
            @inbounds @simd for i in axes(p, 2)
                dst[1] += m[i] * ($p2 * $v3 - $p3 * $v2)
                dst[2] += m[i] * ($p3 * $v1 - $p1 * $v3)
                dst[3] += m[i] * ($p1 * $v2 - $p2 * $v1)
            end
            return dst
        end

        function _angmomtot_stable_unsafe(
            p::AbstractMatrix,
            v::AbstractMatrix,
            m::Number,
            origin::$T1,
            velorigin::$T2,
        )
            return [
                m * sum($p2 * $v3 - $p3 * $v2 for i in axes(p, 2))
                m * sum($p3 * $v1 - $p1 * $v3 for i in axes(p, 2))
                m * sum($p1 * $v2 - $p2 * $v1 for i in axes(p, 2))
            ]
        end

        function _angmomtot_stable_unsafe(
            p::AbstractMatrix,
            v::AbstractMatrix,
            m::AbstractVector,
            origin::$T1,
            velorigin::$T2,
        )
            return [
                sum(m[i] * ($p2 * $v3 - $p3 * $v2) for i in axes(p, 2))
                sum(m[i] * ($p3 * $v1 - $p1 * $v3) for i in axes(p, 2))
                sum(m[i] * ($p1 * $v2 - $p2 * $v1) for i in axes(p, 2))
            ]
        end
    end |> eval
end
