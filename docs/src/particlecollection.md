```@meta
CurrentModule = CosmoParticles
```

# Particle Collections

A set of different particles are represented by subtypes of `AbstractParticleCollection`.
These can have a cosmological redshift attributed to them.
In the context of cosmological simulations, this applies to collections of DM, gas, and star particles.

For an example implementation of a concrete subtype of `AbstractParticleCollection`, see the source codes of
[`ParticleCollection`](@ref) or [`RedshiftParticleCollection`](@ref).

```@docs
AbstractParticleCollection
ParticleCollection
ParticleCollection(::Type{<:AbstractParticles})
ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}
ParticleCollection(p::Particles, ps::Particles...)
RedshiftParticleCollection
RedshiftParticleCollection(::Type{<:AbstractParticles}, z::Real)
RedshiftParticleCollection(z::T, pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles,T<:Real}
RedshiftParticleCollection(z::Real, p::Particles, ps::Particles...)
redshift
CosmoParticles.get_particles
CosmoParticles.show_properties(io::IO, mime::AbstractString, pc::AbstractParticleCollection)
```
