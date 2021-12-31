```@meta
CurrentModule = CosmoParticles
```

# Particles

A set of objects with different properties are represented by subtypes of `AbstractParticles`.
While not absolutely necessary, these kind of objects generally have IDs, positions, and velocity as properties.
In the context of cosmological simulations, this applies to particle and galaxy data.

For an example implementation of a concrete subtype of `AbstractParticles`, see the source code of [`Particles`](@ref).

```@docs
AbstractParticles
Particles
Particles(type)
Particles(p::AllParticles, props=())
CosmoParticles.get_props
CosmoParticles.particle_name
CosmoParticles.particle_number
CosmoParticles.show_properties(io::IO, mime::AbstractString, p::AbstractParticles)
```
