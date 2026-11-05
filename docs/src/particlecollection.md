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
ParticleCollection(particles::AbstractDict{Symbol,AP}) where {AP<:AbstractParticles}
ParticleCollection(::Type{<:AbstractParticles})
ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}
ParticleCollection(p::Particles, ps::Particles...)
RedshiftParticleCollection
RedshiftParticleCollection(z::T, particles::AbstractDict{Symbol,AP}) where {AP<:AbstractParticles,T<:Real}
RedshiftParticleCollection(::Type{<:AbstractParticles}, z::Real)
RedshiftParticleCollection(z::T, pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles,T<:Real}
RedshiftParticleCollection(z::Real, p::Particles, ps::Particles...)
redshift
CosmoParticles.get_particles
CosmoParticles.show_properties(io::IO, mime::AbstractString, pc::AbstractParticleCollection)
```

## Tables

Particle collections implement the [Tables.jl](https://github.com/JuliaData/Tables.jl) interface like particles, with
one row per particle. The first column `ptype` holds the particle type of each particle, and the other columns are
the table columns of all particles as given by [`AllParticles`](@ref). Particles without a property have `missing`
values in its columns.

```jldoctest
julia> using CosmoParticles, Tables

julia> pc = ParticleCollection(
           Particles(:dm, :id => [1, 2], :pos => [0.0 1.0; 0.0 0.0; 0.0 0.0]),
           Particles(:gas, :id => [3], :pos => [2.0; 0.0; 0.0;;], :temp => [1e4]),
       );

julia> Tables.columntable(pc)
(ptype = [:dm, :dm, :gas], id = [1, 2, 3], pos_x = [0.0, 1.0, 2.0], pos_y = [0.0, 0.0, 0.0], pos_z = [0.0, 0.0, 0.0], temp = Union{Missing, Float64}[missing, missing, 10000.0])
```

# All Particles

```@docs
AllParticles
particle_collection
Base.copy!(dst::AbstractParticles, src::AllParticles, props=nothing)
```
