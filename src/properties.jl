"""
    meanprop(p::AbstractParticles, prop::Symbol; massprop=:mass, massweighted=true)

Returns the mean of the specified property over all particles.

The mean is weighted by the mass by default.

# Examples
```jldoctest
julia> p = Particles(:dm, :mass => [1.0, 3.0], :temp => [1.0, 5.0]);

julia> meanprop(p, :temp)
4.0

julia> meanprop(p, :temp; massweighted=false)
3.0
```
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
_meanprop(a::AbstractFillVector) = a[1]
_meanprop(a::AbstractFillVector, ::AbstractVector) = a[1]
_meanprop(a::AbstractFillVector, ::AbstractFillVector) = a[1]
_meanprop(a::AbstractVector) = mean(a)
# NOTE: commented versions had to be removed because the summation over a large number of
# Float32 values with the `dims` keyword calls a mapreduce summation, which can introduce
# large numerical errors. This is not the case when summing without `dims`.
# _meanprop(a::AbstractMatrix) = mean(a; dims=2)
_meanprop(a::AbstractMatrix) = [mean(@view a[i, :]) for i in axes(a, 1)]
_meanprop(a::AbstractVector, m::AbstractVector) = mean(a, weights(ustrip_lazy(m)))
# _meanprop(a::AbstractMatrix, m::AbstractVector) = mean(a, weights(ustrip_lazy(m)); dims=2)
function _meanprop(a::AbstractMatrix, m::AbstractVector)
    w = weights(ustrip_lazy(m))
    return [mean(@view(a[i, :]), w) for i in axes(a, 1)]
end
_meanprop(a::AbstractFillMatrix) = Fill(a[1, 1], size(a, 1))
_meanprop(a::AbstractFillMatrix, ::AbstractVector) = _meanprop(a)
_meanprop(a::AbstractFillMatrix, ::AbstractFillVector) = _meanprop(a)
_meanprop(a::AbstractArray, ::Number) = _meanprop(a)
_meanprop(a::AbstractVector, ::AbstractFillVector) = _meanprop(a)
_meanprop(a::AbstractMatrix, ::AbstractFillVector) = _meanprop(a)


"""
    sumprop(p::AbstractParticles, prop::Symbol)

Returns the sum of the specified property over all particles.

For a matrix property, the sums of its rows are returned as a `Vector`.

# Examples
```jldoctest
julia> p = Particles(:dm, :mass => [1.0, 2.0, 3.0]);

julia> sumprop(p, :mass)
6.0
```
"""
sumprop(p::AbstractParticles, prop::Symbol) = _sumprop(p, p[prop])

_sumprop(p::AbstractParticles, a::Number) = particle_number(p) * a
_sumprop(::AbstractParticles, a::AbstractVecOrMat) = _sumprop(a)
_sumprop(a::AbstractVector) = sum(a)
_sumprop(a::AbstractMatrix) = [sum(@view a[i, :]) for i in axes(a, 1)]


"""
    meanpos(p::AbstractParticles, prop::Symbol=:pos; massprop=:mass, massweighted=true)

Returns the mean of the particles' positions.

The mean is weighted by the mass by default.

# Examples
```jldoctest
julia> p = Particles(:dm, :pos => [0.0 2.0; 0.0 4.0; 0.0 6.0], :mass => [1.0, 1.0]);

julia> meanpos(p)
3-element Vector{Float64}:
 1.0
 2.0
 3.0
```
"""
meanpos(p::AbstractParticles, prop::Symbol=:pos; massprop=:mass, massweighted=true) =
    meanprop(p, prop; massprop, massweighted)

"""
    meanvel(p::AbstractParticles, prop::Symbol=:vel; massprop=:mass, massweighted=true)

Returns the mean of the particles' velocities.

The mean is weighted by the mass by default.

# Examples
```jldoctest
julia> p = Particles(:dm, :vel => [0.0 2.0; 0.0 4.0; 0.0 6.0], :mass => [1.0, 1.0]);

julia> meanvel(p)
3-element Vector{Float64}:
 1.0
 2.0
 3.0
```
"""
meanvel(p::AbstractParticles, prop::Symbol=:vel; massprop=:mass, massweighted=true) =
    meanprop(p, prop; massprop, massweighted)


"""
    angmom(p::AbstractParticles; kwargs...)

Returns the individual angular momenta of the particles.

If `angmomprop` is a `Symbol`, the already computed angular momenta are returned directly if they exist.

# Keyword Arguments
- `origin=nothing`: origin of the coordinate system for the computation of the angular momentum
- `velorigin=nothing`: relative velocity zero-point, oftentimes the velocity of the center of mass
- `posprop=:pos`: property name of the position
- `velprop=:vel`: property name of the velocity
- `massprop=:mass`: property name of the mass
- `angmomprop=nothing`: property name of the angular momentum - the method checks whether it already exists or not. Warning: The method assumes the already saved property to have used the same `origin` and `velorigin`.

# Examples
```jldoctest
julia> pos = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> vel = [0.0 1.0; 1.0 0.0; 0.0 0.0];

julia> p = Particles(:dm, :pos => pos, :vel => vel, :mass => [1.0, 2.0]);

julia> angmom(p)
3×2 Matrix{Float64}:
 0.0   0.0
 0.0   0.0
 1.0  -2.0
```
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
    angmomtot(p::AbstractParticles; kwargs...)

Returns the total angular momentum of the particles.

If `angmomprop` is a `Symbol`, the total angular momentum is computed from the already existing angular momenta.

# Keyword Arguments
- `origin=nothing`: origin of the coordinate system for the computation of the angular momentum
- `velorigin=nothing`: relative velocity zero-point, oftentimes the velocity of the center of mass
- `posprop=:pos`: property name of the position
- `velprop=:vel`: property name of the velocity
- `massprop=:mass`: property name of the mass
- `angmomprop=nothing`: property name of the angular momentum - the method checks whether it already exists or not. Warning: The method assumes the already saved property to have used the same `origin` and `velorigin`.

# Examples
```jldoctest
julia> pos = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> vel = [0.0 1.0; 1.0 0.0; 0.0 0.0];

julia> p = Particles(:dm, :pos => pos, :vel => vel, :mass => [1.0, 2.0]);

julia> angmomtot(p)
3-element Vector{Float64}:
  0.0
  0.0
 -1.0
```
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
    angmomtot_stable(p::AbstractParticles; kwargs...)

Returns the total angular momentum of the particles, like [`angmomtot`](@ref).

Sums with [`Base.sum`](@extref), which adds the terms pairwise and accumulates less rounding error than the running sum
of [`angmomtot`](@ref), but is almost three times slower than it.

See [`angmomtot`](@ref) for the available keyword arguments.

# Examples
```jldoctest
julia> pos = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> vel = [0.0 1.0; 1.0 0.0; 0.0 0.0];

julia> p = Particles(:dm, :pos => pos, :vel => vel, :mass => [1.0, 2.0]);

julia> angmomtot_stable(p)
3-element Vector{Float64}:
  0.0
  0.0
 -1.0
```
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

# Examples
```jldoctest
julia> a = [3.0 0.0 0.0; 4.0 5.0 1.0];

julia> colnorm(a)
3-element Vector{Float64}:
 5.0
 5.0
 1.0

julia> colnorm(a, [3.0, 1.0])
3-element Vector{Float64}:
 3.0
 5.0
 3.0
```
"""
function colnorm(a::AbstractMatrix{T}) where {T}
    dst = Vector{Base.promote_op(float, T)}(undef, size(a, 2))
    return _colnorm_unsafe!(dst, a)
end

function colnorm(a::AbstractMatrix{T}, origin::AbstractVector{T0}) where {T,T0}
    _check_origin_length(a, origin)
    T === T0 || (origin = convert(Vector{T}, origin))
    dst = Vector{Base.promote_op(float, T)}(undef, size(a, 2))
    return _colnorm_unsafe!(dst, a, origin)
end

"""
    colnorm2(a::AbstractMatrix[, origin::AbstractVector])

Returns a new `Vector` with the columnwise squared norms of `a` around `origin`.

# Examples
```jldoctest
julia> a = [3.0 0.0 0.0; 4.0 5.0 1.0];

julia> colnorm2(a)
3-element Vector{Float64}:
 25.0
 25.0
  1.0

julia> colnorm2(a, [0.0, 1.0])
3-element Vector{Float64}:
 18.0
 16.0
  0.0
```
"""
function colnorm2(a::AbstractMatrix{T}) where {T}
    dst = Vector{Base.promote_op(abs2, T)}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a)
end

function colnorm2(a::AbstractMatrix{T}, origin::AbstractVector) where {T}
    _check_origin_length(a, origin)
    dst = Vector{Base.promote_op(abs2, T)}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, _help_convert_origin(T, origin))
end


function _check_origin_length(a::AbstractMatrix, origin::AbstractVector)
    size(a, 1) == length(origin) || throw(
        DimensionMismatch("the origin must have length $(size(a, 1)) like the matrix rows, got length $(length(origin))"),
    )
    return nothing
end

function _check_same_size(a::AbstractMatrix, b::AbstractMatrix)
    size(a) == size(b) || throw(DimensionMismatch("the matrices must have the same size, got $(size(a)) and $(size(b))"))
    return nothing
end

function _colnorm_unsafe! end
function _colnorm2_unsafe! end

for (name, zeroT2) in zip(["_colnorm_unsafe", "_colnorm2_unsafe"], [:(zero(T)^2), :(zero(T))])
    quote
        function $(Symbol(name, "!"))(dst::AbstractVector{T}, a::AbstractMatrix) where {T}
            Base.require_one_based_indexing(dst, a)
            dims = size(a, 1)
            if dims == 2
                return $(Symbol(name, "_2D!"))(dst, a)
            elseif dims == 3
                return $(Symbol(name, "_3D!"))(dst, a)
            else
                return $(Symbol(name, "_ND!"))(dst, a; zeroT2=($zeroT2))
            end
        end

        function $(Symbol(name, "!"))(
            dst::AbstractVector{T},
            a::AbstractMatrix,
            origin::AbstractVector,
        ) where {T}
            Base.require_one_based_indexing(dst, a, origin)
            dims = size(a, 1)
            if dims == 2
                return $(Symbol(name, "_2D!"))(dst, a, origin)
            elseif dims == 3
                return $(Symbol(name, "_3D!"))(dst, a, origin)
            else
                return $(Symbol(name, "_ND!"))(dst, a, origin; zeroT2=($zeroT2))
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

Returns the dot products of the matrix columns as a new vector.

Both input matrices need to have the same dimensions ``d × N``, representing ``N`` ``d``-dimensional vectors.

# Examples
```jldoctest
julia> coldot([1.0 2.0 3.0; 4.0 5.0 6.0], [1.0 0.0 1.0; 0.0 1.0 1.0])
3-element Vector{Float64}:
 1.0
 5.0
 9.0
```
"""
function coldot(a::AbstractMatrix, b::AbstractMatrix)
    _check_same_size(a, b)
    return _coldot_unsafe(a, b)
end

function _coldot_unsafe(a::AbstractMatrix, b::AbstractMatrix)
    T = Base.promote_op((x, y) -> x * y + x * y, eltype(a), eltype(b))
    dst = similar(a, T, size(a, 2))
    return _coldot_unsafe!(dst, a, b)
end

function _coldot_unsafe!(dst::AbstractVector, a::AbstractMatrix, b::AbstractMatrix)
    Base.require_one_based_indexing(dst, a, b)
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

# Examples
```jldoctest
julia> a = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> b = [0.0 0.0; 1.0 0.0; 0.0 1.0];

julia> colcross(a, b)
3×2 Matrix{Float64}:
 0.0  1.0
 0.0  0.0
 1.0  0.0
```
"""
function colcross(a::AbstractMatrix, b::AbstractMatrix)
    _check_same_size(a, b)
    size(a, 1) == 3 || throw(DimensionMismatch("the matrices must have 3 rows, got $(size(a, 1))"))
    return _colcross_unsafe(a, b)
end

function _colcross_unsafe(a::AbstractMatrix, b::AbstractMatrix)
    T = Base.promote_op((x, y) -> x * y - x * y, eltype(a), eltype(b))
    dst = similar(a, T)
    return _colcross_unsafe!(dst, a, b)
end

function _colcross_unsafe!(dst::AbstractMatrix, a::AbstractMatrix, b::AbstractMatrix)
    Base.require_one_based_indexing(dst, a, b)
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

# Examples
```jldoctest
julia> pos = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> vel = [0.0 1.0; 1.0 0.0; 0.0 0.0];

julia> angmom(pos, vel, [1.0, 2.0])
3×2 Matrix{Float64}:
 0.0   0.0
 0.0   0.0
 1.0  -2.0
```
"""
function angmom(
    pos::AbstractMatrix,
    vel::AbstractMatrix,
    mass::AbstractVector;
    origin=nothing,
    velorigin=nothing,
)
    _check_angmom_args(pos, vel, mass, origin, velorigin)
    return _angmom_unsafe(pos, vel, mass, origin, velorigin)
end

function angmom(pos::AbstractMatrix, vel::AbstractMatrix, mass::Number; origin=nothing, velorigin=nothing)
    _check_angmom_args(pos, vel, mass, origin, velorigin)
    return _angmom_unsafe(pos, vel, mass, origin, velorigin)
end


# The angular momentum kernels index the rows of the positions, velocities, and origins as 1 to 3, and the
# particles from 1.
function _check_angmom_args(pos, vel, mass, origin, velorigin)
    size(pos, 1) == 3 || throw(DimensionMismatch("the positions must have 3 rows, got $(size(pos, 1))"))
    size(pos) == size(vel) || throw(
        DimensionMismatch("the positions and velocities must have the same size, got $(size(pos)) and $(size(vel))"),
    )
    mass isa AbstractVector && length(mass) != size(pos, 2) && throw(
        DimensionMismatch("the masses must have length $(size(pos, 2)) like the positions, got length $(length(mass))"),
    )
    Base.require_one_based_indexing(pos, vel)
    mass isa AbstractArray && Base.require_one_based_indexing(mass)
    _check_origin(origin, "origin")
    _check_origin(velorigin, "velorigin")
    return nothing
end

_check_origin(::Nothing, _) = nothing
function _check_origin(origin::AbstractVector, name)
    Base.require_one_based_indexing(origin)
    length(origin) == 3 ||
        throw(DimensionMismatch("`$name` must have length 3, got length $(length(origin))"))
    return nothing
end

# An origin takes the element type of the positions or velocities, without `Missing`.
_help_convert_origin(::Type{T}, origin::AbstractVector) where {T} = convert(Vector{nonmissingtype(T)}, origin)
_help_convert_origin(::Type{Missing}, origin::AbstractVector) = origin
_help_convert_origin(::Type, ::Nothing) = nothing


function _angmom_unsafe(pos::AbstractMatrix, vel::AbstractMatrix, mass, origin, velorigin)
    T = Base.promote_op((x, y, m) -> m * (x * y - x * y), eltype(pos), eltype(vel), _eltype_of(mass))
    dst = similar(pos, T)
    return _angmom_unsafe!(
        dst,
        pos,
        vel,
        mass,
        _help_convert_origin(eltype(pos), origin),
        _help_convert_origin(eltype(vel), velorigin),
    )
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


_eltype_of(m::AbstractArray) = eltype(m)
_eltype_of(m) = typeof(m)



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

# Examples
```jldoctest
julia> pos = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> vel = [0.0 1.0; 1.0 0.0; 0.0 0.0];

julia> angmomtot(pos, vel, [1.0, 2.0])
3-element Vector{Float64}:
  0.0
  0.0
 -1.0
```
"""
function angmomtot(
    pos::AbstractMatrix,
    vel::AbstractMatrix,
    mass::AbstractVector;
    origin=nothing,
    velorigin=nothing,
)
    _check_angmom_args(pos, vel, mass, origin, velorigin)
    return _angmomtot_unsafe(pos, vel, mass, origin, velorigin)
end

function angmomtot(pos::AbstractMatrix, vel::AbstractMatrix, mass::Number; origin=nothing, velorigin=nothing)
    _check_angmom_args(pos, vel, mass, origin, velorigin)
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

# The kernels for a `Fill` mass, like those for a scalar mass, sum without the mass and multiply by it at the end.
function _type_and_zero_of_summed_angmomtot(p::AbstractMatrix, v::AbstractMatrix, m::AbstractFillVector)
    return _type_and_zero_of_summed_angmomtot(p, v, zero(eltype(m)))
end


# The value of a `Fill` mass; an empty `Fill` has none, and its sums are zero.
_fill_value(m::AbstractFillVector) = isempty(m) ? zero(eltype(m)) : m[1]

"""
    angmomtot_stable(pos::AbstractMatrix, vel::AbstractMatrix, mass; origin=nothing, velorigin=nothing)

Returns the total angular momentum based on the positions, velocities, and masses, like [`angmomtot`](@ref).

Sums with [`Base.sum`](@extref), which adds the terms pairwise and accumulates less rounding error than the running sum
of [`angmomtot`](@ref), but is almost three times slower than it.

See [`angmomtot`](@ref) for the available keyword arguments.

# Examples
```jldoctest
julia> using FillArrays

julia> pos = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> vel = [0.0 -1.0; 1.0 0.0; 0.0 0.0];

julia> angmomtot_stable(pos, vel, Fill(2.0, 2))
3-element Vector{Float64}:
 0.0
 0.0
 4.0
```
"""
function angmomtot_stable end

function angmomtot_stable(
    pos::AbstractMatrix,
    vel::AbstractMatrix,
    mass::AbstractVector;
    origin=nothing,
    velorigin=nothing,
)
    _check_angmom_args(pos, vel, mass, origin, velorigin)
    # `sum` throws for no particles, and passing it `init` would make it a running sum.
    size(pos, 2) == 0 && return _angmomtot_unsafe(pos, vel, mass, origin, velorigin)
    return _angmomtot_stable_unsafe(pos, vel, mass, origin, velorigin)
end

function angmomtot_stable(
    pos::AbstractMatrix,
    vel::AbstractMatrix,
    mass::Number;
    origin=nothing,
    velorigin=nothing,
)
    _check_angmom_args(pos, vel, mass, origin, velorigin)
    # `sum` throws for no particles, and passing it `init` would make it a running sum.
    size(pos, 2) == 0 && return _angmomtot_unsafe(pos, vel, mass, origin, velorigin)
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
            m::AbstractFillVector,
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
            mval = _fill_value(m)
            dst .= mval * d1, mval * d2, mval * d3
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
                m * sum(i -> $p2 * $v3 - $p3 * $v2, axes(p, 2))
                m * sum(i -> $p3 * $v1 - $p1 * $v3, axes(p, 2))
                m * sum(i -> $p1 * $v2 - $p2 * $v1, axes(p, 2))
            ]
        end

        function _angmomtot_stable_unsafe(
            p::AbstractMatrix,
            v::AbstractMatrix,
            m::AbstractFillVector,
            origin::$T1,
            velorigin::$T2,
        )
            mval = m[1]
            return [
                mval * sum(i -> $p2 * $v3 - $p3 * $v2, axes(p, 2))
                mval * sum(i -> $p3 * $v1 - $p1 * $v3, axes(p, 2))
                mval * sum(i -> $p1 * $v2 - $p2 * $v1, axes(p, 2))
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
                sum(i -> m[i] * ($p2 * $v3 - $p3 * $v2), axes(p, 2))
                sum(i -> m[i] * ($p3 * $v1 - $p1 * $v3), axes(p, 2))
                sum(i -> m[i] * ($p1 * $v2 - $p2 * $v1), axes(p, 2))
            ]
        end
    end |> eval
end
