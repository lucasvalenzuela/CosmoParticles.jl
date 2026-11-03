"""
    CosmoParticles.matrix_rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})

Returns `vals` rotated by the rotation matrix `rotmat` by regular matrix multiplication.
Works for any dimensions and is optimized for unitful arrays.

# Examples

```jldoctest
julia> CosmoParticles.matrix_rotate([1 0; 0 1; 0 0], [0 -1 0; 1 0 0; 0 0 1])
3×2 Matrix{Int64}:
 0  -1
 1   0
 0   0
```
"""
matrix_rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real}) = rotmat*ᵘvals

# A product keeps the row names of a `NamedRowMatrix` only as its left factor, and the fast unitful product
# applies to the wrapped data, not to the `NamedRowMatrix` itself.
matrix_rotate(vals::NamedRowMatrix{<:Number}, rotmat::AbstractMatrix{<:Real}) =
    NamedRowArray(matrix_rotate(parent(vals), rotmat), rownames(vals))

"""
    CosmoParticles.matrix_rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})

In-place version of [`CosmoParticles.matrix_rotate`](@ref), overwriting `vals`.

# Examples

```jldoctest
julia> vals = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> CosmoParticles.matrix_rotate!(vals, [0 -1 0; 1 0 0; 0 0 1])
3×2 Matrix{Float64}:
 0.0  -1.0
 1.0   0.0
 0.0   0.0
```
"""
function matrix_rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})
    vals .= matrix_rotate(vals, rotmat)
    return vals
end

# Without this special method for lazy arrays, matrix multiplication with unitful arrays leads to bizarre errors
# where only one row or certain diagonal parts of the array or properly rotated.
function matrix_rotate!(vals::ApplyMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})
    matrix_rotate!.(vals.args, (rotmat,))
    return vals
end


"""
    CosmoParticles.matrix_rotate!(dst::AbstractMatrix{<:Number}, vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})

In-place version of [`CosmoParticles.matrix_rotate`](@ref), overwriting `dst`.

# Examples

```jldoctest
julia> vals = [1.0 0.0; 0.0 1.0; 0.0 0.0];

julia> dst = similar(vals);

julia> CosmoParticles.matrix_rotate!(dst, vals, [0 -1 0; 1 0 0; 0 0 1]);

julia> dst
3×2 Matrix{Float64}:
 0.0  -1.0
 1.0   0.0
 0.0   0.0
```
"""
function matrix_rotate!(
    dst::AbstractMatrix{<:Number},
    vals::AbstractMatrix{<:Number},
    rotmat::AbstractMatrix{<:Real},
)
    uproduct!(dst, rotmat, vals)
end


# Explicitly named properties must exist: for particles, or for at least one particle type of a collection.
function _check_has_props(p::AbstractParticles, props)
    for prop in props
        haskey(p, prop) || throw(ArgumentError("the particles do not have the property :$prop"))
    end
    return nothing
end

function _check_has_props(pc::AbstractParticleCollection, props)
    for prop in props
        any(p -> haskey(p, prop), values(pc)) ||
            throw(ArgumentError("no particle type in the collection has the property :$prop"))
    end
    return nothing
end

"""
    rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}[, props=[:pos, :vel]])

Rotates the properties `props` of the particles `p` by the rotation matrix `rotmat`.

Creates a copy of the particles with only new pointers to the rotated properties. This means that mutations to
other properties will affect both the original particles as well as the rotated copied particles.
The properties are given as a single `Symbol` or a vector of `Symbol`s. By default, the positions `:pos` and the
velocities `:vel` are rotated, and those that do not exist are skipped. Properties that are passed must exist.

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [1.0 0.0; 0.0 1.0; 0.0 0.0]);

julia> rotate(p, [0 -1 0; 1 0 0; 0 0 1]).pos
3×2 Matrix{Float64}:
 0.0  -1.0
 1.0   0.0
 0.0   0.0
```
"""
function rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props)
    _check_has_props(p, props)
    return _rotate(p, rotmat, props)
end

rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}) = _rotate(p, rotmat, (:pos, :vel))
rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop::Symbol) = rotate(p, rotmat, (prop,))

# Rotates the properties among `props` that exist.
function _rotate(p::AbstractParticles, rotmat, props)
    p = copy(p)

    Threads.@threads for prop in intersect(keys(p), props) |> collect # collect used for compatibility with threads
        _setprop!(p, prop, matrix_rotate(p[prop], rotmat))
    end

    return p
end

"""
    rotate(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}[, props=[:pos, :vel]])

Rotates the properties `props` of the particles in the collection by the rotation matrix `rotmat`.

Creates a copy of the particle collection with only new pointers to the rotated particles.
The properties are given as a single `Symbol` or a vector of `Symbol`s. By default, the positions `:pos` and the
velocities `:vel` are rotated, and those that do not exist are skipped. Properties that are passed must exist for at
least one particle type.

# Examples

```jldoctest
julia> pc = ParticleCollection(Particles(:dm, :pos => [1.0 0.0; 0.0 1.0; 0.0 0.0]));

julia> rotate(pc, [0 -1 0; 1 0 0; 0 0 1])[:dm].pos
3×2 Matrix{Float64}:
 0.0  -1.0
 1.0   0.0
 0.0   0.0
```
"""
function rotate(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, props)
    _check_has_props(pc, props)
    return _rotate(pc, rotmat, props)
end

rotate(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}) = _rotate(pc, rotmat, (:pos, :vel))
rotate(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, prop::Symbol) = rotate(pc, rotmat, (prop,))

function _rotate(pc::AbstractParticleCollection, rotmat, props)
    pc = copy(pc)

    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        pc[ptype] = _rotate(pc[ptype], rotmat, props)
    end

    return pc
end

"""
    LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}[, props=[:pos, :vel]])
    LinearAlgebra.rotate!(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}[, props=[:pos, :vel]])

In-place version of [`rotate`](@ref).

This function is reexported.

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [1.0 0.0; 0.0 1.0; 0.0 0.0]);

julia> rotate!(p, [0 -1 0; 1 0 0; 0 0 1]);

julia> p.pos
3×2 Matrix{Float64}:
 0.0  -1.0
 1.0   0.0
 0.0   0.0
```
"""
function LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props)
    _check_has_props(p, props)
    return _rotate!(p, rotmat, props)
end

LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}) = _rotate!(p, rotmat, (:pos, :vel))
LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop::Symbol) =
    rotate!(p, rotmat, (prop,))

function _rotate!(p::AbstractParticles, rotmat, props)
    Threads.@threads for prop in intersect(keys(p), props) |> collect # collect used for compatibility with threads
        matrix_rotate!(p[prop], rotmat)
    end

    return p
end

function LinearAlgebra.rotate!(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, props)
    _check_has_props(pc, props)
    return _rotate!(pc, rotmat, props)
end

LinearAlgebra.rotate!(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}) =
    _rotate!(pc, rotmat, (:pos, :vel))
LinearAlgebra.rotate!(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, prop::Symbol) =
    rotate!(pc, rotmat, (prop,))

function _rotate!(pc::AbstractParticleCollection, rotmat, props)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        _rotate!(pc[ptype], rotmat, props)
    end

    return pc
end


"""
    translate(p::AbstractParticles, Δx[, prop::Symbol=:pos])

Translates the property `prop` of the particles `p` by `Δx`.

Creates a copy of the particles with only new pointers to the translated property.
Typically, `Δx` will be an `AbstractVector`, e.g., for positions or velocities.
By default, the positions `:pos` are translated, and skipped if they do not exist. A property that is passed must
exist.

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [1.0 2.0; 0.0 0.0; 0.0 0.0]);

julia> translate(p, [1.0, 2.0, 3.0]).pos
3×2 Matrix{Float64}:
 2.0  3.0
 2.0  2.0
 3.0  3.0
```
"""
function translate(p::AbstractParticles, Δx, prop::Symbol)
    _check_has_props(p, (prop,))
    return _translate(p, Δx, prop)
end

translate(p::AbstractParticles, Δx) = _translate(p, Δx, :pos)

# Translates the property if it exists.
function _translate(p::AbstractParticles, Δx, prop::Symbol)
    p = copy(p)

    if haskey(p, prop)
        _setprop!(p, prop, p[prop] .+ Δx)
    end

    return p
end

"""
    translate(pc::AbstractParticleCollection, Δx[, prop::Symbol=:pos])

Translates the property `prop` of the particles in the collection `pc` by `Δx`.

Creates a copy of the particle collection with only new pointers to the translated properties for the particles.
Typically, `Δx` will be an `AbstractVector`, e.g., for positions or velocities.
By default, the positions `:pos` are translated, and skipped if they do not exist. A property that is passed must
exist for at least one particle type.

# Examples

```jldoctest
julia> pc = ParticleCollection(Particles(:dm, :pos => [1.0 2.0; 0.0 0.0; 0.0 0.0]));

julia> translate(pc, [1.0, 2.0, 3.0])[:dm].pos
3×2 Matrix{Float64}:
 2.0  3.0
 2.0  2.0
 3.0  3.0
```
"""
function translate(pc::AbstractParticleCollection, Δx, prop::Symbol)
    _check_has_props(pc, (prop,))
    return _translate(pc, Δx, prop)
end

translate(pc::AbstractParticleCollection, Δx) = _translate(pc, Δx, :pos)

function _translate(pc::AbstractParticleCollection, Δx, prop::Symbol)
    pc = copy(pc)

    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        pc[ptype] = _translate(pc[ptype], Δx, prop)
    end

    return pc
end

"""
    translate!(p::AbstractParticles, Δx[, prop::Symbol=:pos])
    translate!(pc::AbstractParticleCollection, Δx[, prop::Symbol=:pos])

In-place version of [`translate`](@ref).

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [1.0 2.0; 0.0 0.0; 0.0 0.0]);

julia> translate!(p, [1.0, 2.0, 3.0]);

julia> p.pos
3×2 Matrix{Float64}:
 2.0  3.0
 2.0  2.0
 3.0  3.0
```
"""
function translate!(p::AbstractParticles, Δx, prop::Symbol)
    _check_has_props(p, (prop,))
    return _translate!(p, Δx, prop)
end

translate!(p::AbstractParticles, Δx) = _translate!(p, Δx, :pos)

# Translates the property if it exists.
function _translate!(p::AbstractParticles, Δx, prop::Symbol)
    if haskey(p, prop)
        p[prop] .+= Δx
    end

    return p
end

function translate!(pc::AbstractParticleCollection, Δx, prop::Symbol)
    _check_has_props(pc, (prop,))
    return _translate!(pc, Δx, prop)
end

translate!(pc::AbstractParticleCollection, Δx) = _translate!(pc, Δx, :pos)

function _translate!(pc::AbstractParticleCollection, Δx, prop::Symbol)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        _translate!(pc[ptype], Δx, prop)
    end

    return pc
end


"""
    CosmoParticles._translate_periodic(x::Number, x₀::Number, period::Number, period_half::Number=1//2*period)

Returns the coordinate `x` shifted by the `period` if reference coordinate `x₀` is on the other side.

For periodic simulation boxes, the period corresponds to the boxlength.
"""
@inline function _translate_periodic(
    x::Number,
    x₀::Number,
    period::Number,
    period_half::Number=1 // 2 * period,
)
    if x - x₀ > period_half
        return x - period
    elseif x₀ - x > period_half
        return x + period
    end
    return x
end

# The center and the period are converted to the element type of the coordinates, so that every branch of
# `_translate_periodic` returns that type and the translated coordinates keep it.
function _periodic_args(vals::AbstractArray, centerpos, period)
    T = nonmissingtype(eltype(vals))
    centerpos, period = _convert_elements(T, centerpos), _convert_elements(T, period)
    return centerpos, period, period ./ 2
end

_convert_elements(::Type{T}, x::Number) where {T} = convert(T, x)
_convert_elements(::Type{T}, x::AbstractArray) where {T} = convert.(T, x)


"""
    translate_periodic(p::AbstractParticles, centerpos, period[, prop::Symbol=:pos])
    translate_periodic(pc::AbstractParticleCollection, centerpos, period[, prop::Symbol=:pos])

Translates the specified property of the particles or collection by the `period` such that all particles are the closest to `centerpos` as possible.

Creates a copy of the particles or the collection with only new pointers to the translated property.
Typically, `centerpos` will be an `AbstractVector`, e.g., for positions, and `period` can be a `Number` or
`AbstractVector`, depending on if the period is the same for all dimensions or if they are different.
By default, the positions `:pos` are translated, and skipped if they do not exist. A property that is passed must
exist for the particles, or for at least one particle type of a collection.

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [0.5 9.5; 0.0 0.0; 0.0 0.0]);

julia> translate_periodic(p, [0.0, 0.0, 0.0], 10.0).pos
3×2 Matrix{Float64}:
 0.5  -0.5
 0.0   0.0
 0.0   0.0
```
"""
function translate_periodic(p::AbstractParticles, centerpos, period, prop::Symbol)
    _check_has_props(p, (prop,))
    return _periodic(p, centerpos, period, prop)
end

translate_periodic(p::AbstractParticles, centerpos, period) = _periodic(p, centerpos, period, :pos)

# Translates the property if it exists.
function _periodic(p::AbstractParticles, centerpos, period, prop::Symbol)
    p = copy(p)

    if haskey(p, prop)
        c, P, Phalf = _periodic_args(p[prop], centerpos, period)
        _setprop!(p, prop, _translate_periodic.(p[prop], c, P, Phalf))
    end

    return p
end

function translate_periodic(pc::AbstractParticleCollection, centerpos, period, prop::Symbol)
    _check_has_props(pc, (prop,))
    return _periodic(pc, centerpos, period, prop)
end

translate_periodic(pc::AbstractParticleCollection, centerpos, period) = _periodic(pc, centerpos, period, :pos)

function _periodic(pc::AbstractParticleCollection, centerpos, period, prop::Symbol)
    pc = copy(pc)

    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        pc[ptype] = _periodic(pc[ptype], centerpos, period, prop)
    end

    return pc
end

"""
    translate_periodic!(p::AbstractParticles, centerpos, period[, prop::Symbol=:pos])
    translate_periodic!(pc::AbstractParticleCollection, centerpos, period[, prop::Symbol=:pos])

In-place version of [`translate_periodic`](@ref).

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [0.5 9.5; 0.0 0.0; 0.0 0.0]);

julia> translate_periodic!(p, [0.0, 0.0, 0.0], 10.0);

julia> p.pos
3×2 Matrix{Float64}:
 0.5  -0.5
 0.0   0.0
 0.0   0.0
```
"""
function translate_periodic!(p::AbstractParticles, centerpos, period, prop::Symbol)
    _check_has_props(p, (prop,))
    return _periodic!(p, centerpos, period, prop)
end

translate_periodic!(p::AbstractParticles, centerpos, period) = _periodic!(p, centerpos, period, :pos)

# Translates the property if it exists.
function _periodic!(p::AbstractParticles, centerpos, period, prop::Symbol)
    if haskey(p, prop)
        c, P, Phalf = _periodic_args(p[prop], centerpos, period)
        p[prop] .= _translate_periodic.(p[prop], c, P, Phalf)
    end

    return p
end

function translate_periodic!(pc::AbstractParticleCollection, centerpos, period, prop::Symbol)
    _check_has_props(pc, (prop,))
    return _periodic!(pc, centerpos, period, prop)
end

translate_periodic!(pc::AbstractParticleCollection, centerpos, period) = _periodic!(pc, centerpos, period, :pos)

function _periodic!(pc::AbstractParticleCollection, centerpos, period, prop::Symbol)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        _periodic!(pc[ptype], centerpos, period, prop)
    end

    return pc
end


"""
    translate_periodic_to_center(p::AbstractParticles, centerpos, period[, prop::Symbol=:pos])
    translate_periodic_to_center(pc::AbstractParticleCollection, centerpos, period[, prop::Symbol=:pos])

Translates the specified property of the particles or collection to place `centerpos` at the origin.

The particles are translated periodically to place them closest to `centerpos`.

Creates a copy of the particles or the collection with only new pointers to the translated property.
Typically, `centerpos` will be an `AbstractVector`, e.g., for positions, and `period` can be a `Number` or
`AbstractVector`, depending on if the period is the same for all dimensions or if they are different.
By default, the positions `:pos` are translated, and skipped if they do not exist. A property that is passed must
exist for the particles, or for at least one particle type of a collection.

Also see [`translate_periodic`](@ref).

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [10.5 19.5; 0.0 0.0; 0.0 0.0]);

julia> translate_periodic_to_center(p, [10.0, 0.0, 0.0], 10.0).pos
3×2 Matrix{Float64}:
 0.5  -0.5
 0.0   0.0
 0.0   0.0
```
"""
function translate_periodic_to_center(p::AbstractParticles, centerpos, period, prop::Symbol)
    _check_has_props(p, (prop,))
    return _periodic_to_center(p, centerpos, period, prop)
end

translate_periodic_to_center(p::AbstractParticles, centerpos, period) = _periodic_to_center(p, centerpos, period, :pos)

# Translates the property if it exists.
function _periodic_to_center(p::AbstractParticles, centerpos, period, prop::Symbol)
    p = copy(p)

    if haskey(p, prop)
        c, P, Phalf = _periodic_args(p[prop], centerpos, period)
        _setprop!(p, prop, _translate_periodic.(p[prop], c, P, Phalf) .- c)
    end

    return p
end

function translate_periodic_to_center(pc::AbstractParticleCollection, centerpos, period, prop::Symbol)
    _check_has_props(pc, (prop,))
    return _periodic_to_center(pc, centerpos, period, prop)
end

translate_periodic_to_center(pc::AbstractParticleCollection, centerpos, period) =
    _periodic_to_center(pc, centerpos, period, :pos)

function _periodic_to_center(pc::AbstractParticleCollection, centerpos, period, prop::Symbol)
    pc = copy(pc)

    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        pc[ptype] = _periodic_to_center(pc[ptype], centerpos, period, prop)
    end

    return pc
end

"""
    translate_periodic_to_center!(p::AbstractParticles, centerpos, period[, prop::Symbol=:pos])
    translate_periodic_to_center!(pc::AbstractParticleCollection, centerpos, period[, prop::Symbol=:pos])

In-place version of [`translate_periodic_to_center`](@ref).

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [10.5 19.5; 0.0 0.0; 0.0 0.0]);

julia> translate_periodic_to_center!(p, [10.0, 0.0, 0.0], 10.0);

julia> p.pos
3×2 Matrix{Float64}:
 0.5  -0.5
 0.0   0.0
 0.0   0.0
```
"""
function translate_periodic_to_center!(p::AbstractParticles, centerpos, period, prop::Symbol)
    _check_has_props(p, (prop,))
    return _periodic_to_center!(p, centerpos, period, prop)
end

translate_periodic_to_center!(p::AbstractParticles, centerpos, period) =
    _periodic_to_center!(p, centerpos, period, :pos)

# Translates the property if it exists.
function _periodic_to_center!(p::AbstractParticles, centerpos, period, prop::Symbol)
    if haskey(p, prop)
        c, P, Phalf = _periodic_args(p[prop], centerpos, period)
        p[prop] .= _translate_periodic.(p[prop], c, P, Phalf) .- c
    end

    return p
end

function translate_periodic_to_center!(pc::AbstractParticleCollection, centerpos, period, prop::Symbol)
    _check_has_props(pc, (prop,))
    return _periodic_to_center!(pc, centerpos, period, prop)
end

translate_periodic_to_center!(pc::AbstractParticleCollection, centerpos, period) =
    _periodic_to_center!(pc, centerpos, period, :pos)

function _periodic_to_center!(pc::AbstractParticleCollection, centerpos, period, prop::Symbol)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        _periodic_to_center!(pc[ptype], centerpos, period, prop)
    end

    return pc
end


"""
    CosmoParticles.factor_to_comoving(z::Real, n::Real)

Returns ``(1 + z)^n``.

# Examples

```jldoctest
julia> CosmoParticles.factor_to_comoving(1, 2)
4.0
```
"""
function factor_to_comoving(z::Real, n::Real)
    float(1 + z)^n
end

"""
    CosmoParticles.factor_to_physical(z::Real, n::Real)

Returns ``(1 + z)^{-n}``.

# Examples

```jldoctest
julia> CosmoParticles.factor_to_physical(1, 2)
0.25
```
"""
function factor_to_physical(z::Real, n::Real)
    1 / float(1 + z)^n
end


"""
    to_comoving(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
    to_comoving(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])

Create new particles or collection with particle properties converted from physical to comoving.

The properties and the positional exponent of their units (e.g., 1 for positions, 3 for volumes, and -3 for
densities) are passed as the keyword argument `propexp` as a vector of tuples. By default, the positions `:pos`
and the velocities `:vel` are converted with an exponent of 1. Properties in `propexp` that do not exist are
skipped.
The redshift `z` is obtained from [`redshift`](@ref) for the particle collection.

The properties are multiplied by ``(1 + z)^n`` where ``n`` is the positional exponent.
The particles do not record whether their properties are comoving or physical, so the properties are
converted on every call and must be physical beforehand.

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [1.0 2.0; 0.0 0.0; 0.0 0.0]);

julia> to_comoving(p, 1).pos
3×2 Matrix{Float64}:
 2.0  4.0
 0.0  0.0
 0.0  0.0
```
"""
function to_comoving end

"""
    to_comoving!(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
    to_comoving!(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])

Convert particle properties from physical to comoving in-place.

In-place version of [`to_comoving`](@ref).

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [1.0 2.0; 0.0 0.0; 0.0 0.0]);

julia> to_comoving!(p, 1);

julia> p.pos
3×2 Matrix{Float64}:
 2.0  4.0
 0.0  0.0
 0.0  0.0
```
"""
function to_comoving! end

"""
    to_physical(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
    to_physical(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])

Create new particles or collection with particle properties converted from comoving to physical.

The properties and the positional exponent of their units (e.g., 1 for positions, 3 for volumes, and -3 for
densities) are passed as the keyword argument `propexp` as a vector of tuples. By default, the positions `:pos`
and the velocities `:vel` are converted with an exponent of 1. Properties in `propexp` that do not exist are
skipped.
The redshift `z` is obtained from [`redshift`](@ref) for the particle collection.

The properties are multiplied by ``(1 + z)^{-n}`` where ``n`` is the positional exponent.
The particles do not record whether their properties are comoving or physical, so the properties are
converted on every call and must be comoving beforehand.

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [2.0 4.0; 0.0 0.0; 0.0 0.0]);

julia> to_physical(p, 1).pos
3×2 Matrix{Float64}:
 1.0  2.0
 0.0  0.0
 0.0  0.0
```
"""
function to_physical end

"""
    to_physical!(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
    to_physical!(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])

Convert particle properties from comoving to physical in-place.

In-place version of [`to_physical`](@ref).

# Examples

```jldoctest
julia> p = Particles(:dm, :pos => [2.0 4.0; 0.0 0.0; 0.0 0.0]);

julia> to_physical!(p, 1);

julia> p.pos
3×2 Matrix{Float64}:
 1.0  2.0
 0.0  0.0
 0.0  0.0
```
"""
function to_physical! end

for (name, factorexpr) in zip(["to_comoving", "to_physical"], [:(float(1 + z)), :(1 / (1 + z))])
    quote
        function $(Symbol(name))(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
            p = copy(p)

            factor = $factorexpr
            Threads.@threads for (prop, n) in propexp
                if haskey(p, prop)
                    _setprop!(p, prop, product_preserve_type(p[prop], factor^n))
                end
            end

            return p
        end

        function $(Symbol(name))(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])
            pc = copy(pc)

            Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
                pc[ptype] = $(Symbol(name))(pc[ptype], redshift(pc); propexp)
            end

            return pc
        end

        function $(Symbol(name, "!"))(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
            factor = $factorexpr
            Threads.@threads for (prop, n) in propexp
                if haskey(p, prop)
                    product_preserve_type!(p[prop], factor^n)
                end
            end

            return p
        end

        function $(Symbol(name, "!"))(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])
            Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
                $(Symbol(name, "!"))(pc[ptype], redshift(pc); propexp)
            end

            return pc
        end
    end |> eval
end
