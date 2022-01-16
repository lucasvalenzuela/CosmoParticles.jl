"""
    CosmoParticles.matrix_rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})

Returns `vals` rotated by the rotation matrix `rotmat` by regular matrix multiplication.
Works for any dimensions and is optimized for unitful arrays.

This is not exported.
"""
matrix_rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real}) = rotmat *ᵘ vals

"""
    CosmoParticles.matrix_rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})

In-place version of [`CosmoParticles.matrix_rotate`](@ref).

This is not exported.
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
    rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop)
    rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props=[:pos, :vel])

Rotates the specified properties `props` of the particles `p` by the rotation matrix `rotmat`.

Creates a copy of the particles with only new pointers to the rotated properties
(by default position and velocity). This means that mutations to other properties will
affect both the original particles as well as the rotated copied particles.
The properties should be given as a single `Symbol`, or a vector of `Symbol`s. Only
the existing quantities will be rotated (e.g., this method does not throw an error if
velocities are not given for the particles).
"""
function rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props=[:pos, :vel])
    p = copy(p)

    # only rotate existing quantities
    for prop in intersect(keys(p), props)
        p[prop] = matrix_rotate(p[prop], rotmat)
    end

    return p
end

rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop::Symbol) = rotate(p, rotmat, (prop,))

"""
    rotate(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, prop)
    rotate(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, props=[:pos, :vel])

Rotates the specified properties `props` of the particles in the collection by the rotation matrix `rotmat`.

Creates a copy of the particle collection with only new pointers to the rotated particles
(by default position and velocity). The properties should be given as a single `Symbol`, or a vector of
`Symbol`s.
"""
function rotate(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, props=[:pos, :vel])
    pc = copy(pc)

    for ptype in keys(pc)
        pc[ptype] = rotate(pc[ptype], rotmat, props)
    end

    return pc
end

function rotate(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, prop::Symbol)
    rotate(pc, rotmat, (prop,))
end

"""
    LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props=[:pos, :vel])
    LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop)
    LinearAlgebra.rotate!(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, props=[:pos, :vel])
    LinearAlgebra.rotate!(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, prop)

In-place version of [`rotate`](@ref).

This function is reexported.
"""
function LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props=[:pos, :vel])
    # only rotate existing quantities
    for prop in intersect(keys(p), props)
        matrix_rotate!(p[prop], rotmat)
    end

    return p
end

LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop::Symbol) =
    rotate!(p, rotmat, (prop,))

function LinearAlgebra.rotate!(
    pc::AbstractParticleCollection,
    rotmat::AbstractMatrix{<:Real},
    props=[:pos, :vel],
)
    for ptype in keys(pc)
        rotate!(pc[ptype], rotmat, props)
    end

    return pc
end

LinearAlgebra.rotate!(pc::AbstractParticleCollection, rotmat::AbstractMatrix{<:Real}, prop::Symbol) =
    rotate!(pc, rotmat, (prop,))


"""
    translate(p::AbstractParticles, Δx, prop::Symbol=:pos)

Translates the specified property of the particles `p` by `Δx`.

Creates a copy of the particles with only new pointers to the translated property.
Typically, `Δx` will be an `AbstractVector`, e.g., for positions or velocities.
"""
function translate(p::AbstractParticles, Δx::AbstractVector{<:Number}, prop::Symbol=:pos)
    p = copy(p)

    if haskey(p, prop)
        p[prop] = p[prop] .+ Δx
    end

    return p
end

"""
    translate(pc::AbstractParticleCollection, Δx::AbstractVector{<:Number}, prop::Symbol=:pos)

Translates the specified property of the particles in the collection `pc` by `Δx`.

Creates a copy of the particle collection with only new pointers to the translated properties for the particles.
Typically, `Δx` will be an `AbstractVector`, e.g., for positions or velocities.
"""
function translate(pc::AbstractParticleCollection, Δx::AbstractVector{<:Number}, prop::Symbol=:pos)
    pc = copy(pc)

    for ptype in keys(pc)
        pc[ptype] = translate(pc[ptype], Δx, prop)
    end

    return pc
end

"""
    translate!(p::AbstractParticles, Δx, prop::Symbol=:pos)
    translate!(pc::AbstractParticleCollection, Δx, prop::Symbol=:pos)

In-place version of [`translate`](@ref).
"""
function translate!(p::AbstractParticles, Δx::AbstractVector{<:Number}, prop::Symbol=:pos)
    if haskey(p, prop)
        p[prop] .+= Δx
    end

    return p
end

function translate!(pc::AbstractParticleCollection, Δx::AbstractVector{<:Number}, prop::Symbol=:pos)
    for ptype in keys(pc)
        translate!(pc[ptype], Δx, prop)
    end

    return pc
end


"""
    to_comoving(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
    to_comoving(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])

Create new particles or collection with particle properties converted from physical to comoving.

The properties and the positional exponent of their units (e.g., 1 for positions, 3 for volumes, and -3 for
densities) are passed as the keyword argument `propexp` as a vector of tuples.
The redshift `z` is obtained from [`redshift`](@ref) for the particle collection.

The properties are multiplied by ``(1 + z)^n`` where ``n`` is the positional exponent.
"""
function to_comoving end

"""
    to_comoving!(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
    to_comoving!(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])

Convert particle properties from physical to comoving in-place.

In-place version of [`to_comoving`](@ref).
"""
function to_comoving! end

"""
    to_physical(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
    to_physical(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])

Create new particles or collection with particle properties converted from comoving to physical.

The properties and the positional exponent of their units (e.g., 1 for positions, 3 for volumes, and -3 for
densities) are passed as the keyword argument `propexp` as a vector of tuples.
The redshift `z` is obtained from [`redshift`](@ref) for the particle collection.

The properties are multiplied by ``(1 + z)^{-n}`` where ``n`` is the positional exponent.
"""
function to_physical end

"""
    to_physical!(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
    to_physical!(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])

Convert particle properties from comoving to physical in-place.

In-place version of [`to_physical`](@ref).
"""
function to_physical! end

for (name, factorexpr) in zip(["to_comoving", "to_physical"], [:(1 + z), :(1 / (1 + z))])
    quote
        function $(Symbol(name))(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
            p = copy(p)

            factor = $factorexpr
            for (prop, n) in propexp
                if haskey(p, prop)
                    p[prop] = product_preserve_type(p[prop], factor^n)
                end
            end

            return p
        end

        function $(Symbol(name))(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])
            pc = copy(pc)

            for ptype in keys(pc)
                pc[ptype] = $(Symbol(name))(pc[ptype], redshift(pc); propexp)
            end

            return pc
        end

        function $(Symbol(name, "!"))(p::AbstractParticles, z::Real; propexp=[(:pos, 1), (:vel, 1)])
            factor = $factorexpr
            for (prop, n) in propexp
                if haskey(p, prop)
                    product_preserve_type!(p[prop], factor^n)
                end
            end

            return p
        end

        function $(Symbol(name, "!"))(pc::AbstractParticleCollection; propexp=[(:pos, 1), (:vel, 1)])
            for ptype in keys(pc)
                $(Symbol(name, "!"))(pc[ptype], redshift(pc); propexp)
            end

            return pc
        end
    end |> eval
end
