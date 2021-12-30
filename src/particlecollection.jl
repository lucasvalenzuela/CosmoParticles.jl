"""
    struct ParticleCollection{AP<:AbstractParticles} <: AbstractParticleCollection{AP}
        particles::Dict{Symbol,AP}
    end

Basic particle collection at redshift 0.
"""
struct ParticleCollection{AP<:AbstractParticles} <: AbstractParticleCollection{AP}
    particles::Dict{Symbol,AP}
end

"""
    ParticleCollection(::Type{AP}=Particles) where {AP<:AbstractParticles}

Create an empty particle collection of the passed type.

By default, the particle collection is created for the [`Particles`] type.
"""
ParticleCollection(::Type{AP}=Particles) where {AP<:AbstractParticles} =
    ParticleCollection{AP}(Dict{Symbol,AP}())

"""
    ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}

Create a particle collection from pairs of Symbols and particles.

# Example
```julia
ParticleCollection(:dm => Particles(:dm), :gas => Particles(:gas))
```
"""
function ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}
    ParticleCollection{AP}(Dict{Symbol,AP}(pairs...))
end

"""
    ParticleCollection(p::Particles...)

Create a particle collection from multiple [`Particles`](@ref).

The `Dict` keys are taken from the [`Particles`](@ref) type.

# Example
```julia
ParticleCollection(Particles(:dm), Particles(:gas))
```
"""
function ParticleCollection(p::Particles, ps::Particles...)
    particles = Dict(p.type => p)
    for particle in ps
        particles[particle.type] = particle
    end
    return ParticleCollection(particles)
end

Base.copy(pc::ParticleCollection) = ParticleCollection(copy(pc.particles))
Base.empty(pc::ParticleCollection) = ParticleCollection(empty(pc.particles))
Base.:(==)(pc1::ParticleCollection, pc2::ParticleCollection) = pc1.particles == pc2.particles

function Base.show(io::IO, mime::MIME"text/plain", pc::ParticleCollection)
    printstyled(io, "ParticleCollection"; bold=true)
    println(io)
    show_properties(io, mime, pc)
end


"""
    struct RedshiftParticleCollection{AP<:AbstractParticles,T<:Real} <: AbstractParticleCollection{AP}
        z::T
        particles::Dict{Symbol,AP}
    end

Particle collection at a given redshift.
"""
struct RedshiftParticleCollection{AP<:AbstractParticles,T<:Real} <: AbstractParticleCollection{AP}
    z::T
    particles::Dict{Symbol,AP}
end

"""
    RedshiftParticleCollection([::Type{AP}=Particles,] z::T) where {AP<:AbstractParticles,T<:Real}

Create an empty particle collection of the passed type at redshift `z`.

By default, the particle collection is created for the [`Particles`] type.
"""
function RedshiftParticleCollection(::Type{AP}, z::T) where {AP<:AbstractParticles,T<:Real}
    return RedshiftParticleCollection{AP,T}(z, Dict{Symbol,AP}())
end
function RedshiftParticleCollection(z::T) where {T<:Real}
    return RedshiftParticleCollection(Particles, z)
end

"""
    RedshiftParticleCollection(z::T, pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles,T<:Real}

Create a particle collection at redshift `z` from pairs of Symbols and particles.

# Example
```julia
RedshiftParticleCollection(0.5, :dm => Particles(:dm), :gas => Particles(:gas))
```
"""
function RedshiftParticleCollection(z::T, pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles,T<:Real}
    RedshiftParticleCollection{AP,T}(z, Dict{Symbol,AP}(pairs...))
end

"""
    RedshiftParticleCollection(z::Real, p::Particles...)

Create a particle collection at redshift `z` from multiple [`Particles`](@ref).

The `Dict` keys are taken from the [`Particles`](@ref) type.

# Example
```julia
RedshiftParticleCollection(0.5, Particles(:dm), Particles(:gas))
```
"""
function RedshiftParticleCollection(z::Real, p::Particles, ps::Particles...)
    particles = Dict(p.type => p)
    for particle in ps
        particles[particle.type] = particle
    end
    return RedshiftParticleCollection(z, particles)
end

Base.copy(pc::RedshiftParticleCollection) = RedshiftParticleCollection(copy(pc.z), copy(pc.particles))
Base.empty(pc::RedshiftParticleCollection) = RedshiftParticleCollection(pc.z, empty(pc.particles))
function Base.:(==)(pc1::RedshiftParticleCollection, pc2::RedshiftParticleCollection)
    return pc1.z == pc2.z && pc1.particles == pc2.particles
end

function Base.propertynames(pc::RedshiftParticleCollection; private=false)
    if private
        return [fieldnames(RedshiftParticleCollection) |> collect; :all; keys(pc) |> collect]
    else
        return [:z; :all; keys(pc) |> collect]
    end
end

redshift(pc::RedshiftParticleCollection) = pc.z

function Base.show(io::IO, mime::MIME"text/plain", pc::RedshiftParticleCollection)
    printstyled(io, "ParticleCollection"; bold=true)
    println(io, " at z = $(round(redshift(pc); sigdigits=4))")
    show_properties(io, mime, pc)
end


function Base.:(==)(pc::ParticleCollection, rpc::RedshiftParticleCollection)
    return rpc.z == 0 && pc.particles == rpc.particles
end
Base.:(==)(rpc::RedshiftParticleCollection, pc::ParticleCollection) = ==(pc, rpc)
