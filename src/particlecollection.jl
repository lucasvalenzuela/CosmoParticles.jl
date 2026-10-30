"""
    struct ParticleCollection{AP<:AbstractParticles} <: AbstractParticleCollection{AP}
        particles::OrderedDict{Symbol,AP}
    end

Basic particle collection at redshift 0.
"""
struct ParticleCollection{AP<:AbstractParticles} <: AbstractParticleCollection{AP}
    particles::OrderedDict{Symbol,AP}

    function ParticleCollection{AP}(particles::AbstractDict) where {AP<:AbstractParticles}
        return new{AP}(_ordered_dict(Symbol, AP, particles))
    end
end

"""
    ParticleCollection(particles::AbstractDict{Symbol,AP}) where {AP<:AbstractParticles}

Create a particle collection from a dictionary of particles.

An `OrderedDict{Symbol,AP}` is stored as it is, and other dictionaries are copied into a new `OrderedDict`.
"""
function ParticleCollection(particles::AbstractDict{Symbol,AP}) where {AP<:AbstractParticles}
    return ParticleCollection{AP}(particles)
end

"""
    ParticleCollection(::Type{AP}=Particles) where {AP<:AbstractParticles}

Create an empty particle collection of the passed type.

By default, the particle collection is created for the [`Particles`] type.
"""
ParticleCollection(::Type{AP}=Particles) where {AP<:AbstractParticles} =
    ParticleCollection{AP}(OrderedDict{Symbol,AP}())

"""
    ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}

Create a particle collection from pairs of Symbols and particles.

# Example
```julia
ParticleCollection(:dm => Particles(:dm), :gas => Particles(:gas))
```
"""
function ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}
    ParticleCollection{AP}(OrderedDict{Symbol,AP}(pairs...))
end

"""
    ParticleCollection(p::Particles...)

Create a particle collection from multiple [`Particles`](@ref).

The keys are taken from the [`Particles`](@ref) type.

# Example
```julia
ParticleCollection(Particles(:dm), Particles(:gas))
```
"""
function ParticleCollection(p::Particles, ps::Particles...)
    particles = OrderedDict(p.type => p)
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
        particles::OrderedDict{Symbol,AP}
    end

Particle collection at a given redshift.
"""
struct RedshiftParticleCollection{AP<:AbstractParticles,T<:Real} <: AbstractParticleCollection{AP}
    z::T
    particles::OrderedDict{Symbol,AP}

    function RedshiftParticleCollection{AP,T}(z, particles::AbstractDict) where {AP<:AbstractParticles,T<:Real}
        return new{AP,T}(z, _ordered_dict(Symbol, AP, particles))
    end
end

"""
    RedshiftParticleCollection(z::T, particles::AbstractDict{Symbol,AP}) where {AP<:AbstractParticles,T<:Real}

Create a particle collection at redshift `z` from a dictionary of particles.

An `OrderedDict{Symbol,AP}` is stored as it is, and other dictionaries are copied into a new `OrderedDict`.
"""
function RedshiftParticleCollection(z::T, particles::AbstractDict{Symbol,AP}) where {AP<:AbstractParticles,T<:Real}
    return RedshiftParticleCollection{AP,T}(z, particles)
end

"""
    RedshiftParticleCollection([::Type{AP}=Particles,] z::T) where {AP<:AbstractParticles,T<:Real}

Create an empty particle collection of the passed type at redshift `z`.

By default, the particle collection is created for the [`Particles`] type.
"""
function RedshiftParticleCollection(::Type{AP}, z::T) where {AP<:AbstractParticles,T<:Real}
    return RedshiftParticleCollection{AP,T}(z, OrderedDict{Symbol,AP}())
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
    RedshiftParticleCollection{AP,T}(z, OrderedDict{Symbol,AP}(pairs...))
end

"""
    RedshiftParticleCollection(z::Real, p::Particles...)

Create a particle collection at redshift `z` from multiple [`Particles`](@ref).

The keys are taken from the [`Particles`](@ref) type.

# Example
```julia
RedshiftParticleCollection(0.5, Particles(:dm), Particles(:gas))
```
"""
function RedshiftParticleCollection(z::Real, p::Particles, ps::Particles...)
    particles = OrderedDict(p.type => p)
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

function Base.propertynames(pc::RedshiftParticleCollection, private::Bool=false)
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
