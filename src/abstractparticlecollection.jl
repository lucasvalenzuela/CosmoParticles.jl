"""
    abstract type AbstractParticleCollection{AP} end

Abstract supertype for storing data of multiple particle types.

The particles are expected to be stored in a `Dict{Symbol,<:AbstractParticles}`.

The inbuilt functionality of `AbstractParticleCollection` includes accessing and setting properties
via the following syntax:
```julia
P<:AbstractParticles
pc::<:AbstractParticleCollection{P}
p::<:P
pc.dm = p
pc[:dm] === p
```

If the struct has additional fields, these can also be accessed by `p.field`, but not by
`p[:field]`. The latter syntax can only be used to access the particles.

# Methods
The methods `Base.keys`, `Base.values`, `Base.haskey`, `Base.empty`, `Base.empty!`, and `Base.isempty`
are forwarded to the property `Dict`.

Concrete types of `AbstractParticleCollection` should have the following methods implemented
(also see the implementation of [`ParticleCollection`](@ref)):
- `Base.copy`: returns new object containing a copy of the `Dict`
- `Base.empty`: returns an identical object, but with an empty `Dict`
- `Base.:(==)`
- [`redshift`](@ref): implement this if the redshift of the particle collection is non-zero
- `Base.propertynames`: implement this if there are additional struct fields
- `Base.show(io, mime, pc)`
"""
abstract type AbstractParticleCollection{AP} end

"""
    get_particles(pc::AbstractParticleCollection)

Return the particles `Dict` belonging to the particle collection.

This returns `pc.particles` by default if not overridden.

This is not exported.
"""
get_particles(pc::AbstractParticleCollection) = pc.particles

function Base.getproperty(pc::APC, sym::Symbol) where {APC<:AbstractParticleCollection}
    if sym in fieldnames(APC)
        return getfield(pc, sym)
    else
        return getindex(pc, sym)
    end
end

function Base.setproperty!(pc::APC, sym::Symbol, val) where {APC<:AbstractParticleCollection}
    if sym in fieldnames(APC)
        setfield!(pc, sym, val)
    else
        setindex!(pc, val, sym)
    end
end

Base.getindex(pc::AbstractParticleCollection, sym::Symbol) = getindex(get_particles(pc), sym)
Base.setindex!(pc::AbstractParticleCollection, val, sym::Symbol) =
    (setindex!(get_particles(pc), val, sym); pc)
Base.keys(pc::AbstractParticleCollection) = keys(get_particles(pc))
Base.haskey(pc::AbstractParticleCollection, key) = key in keys(pc)
Base.values(pc::AbstractParticleCollection) = values(get_particles(pc))
Base.propertynames(pc::AbstractParticleCollection) = keys(pc) |> collect
Base.empty!(pc::AbstractParticleCollection) = (empty!(get_particles(pc)); pc)
Base.isempty(pc::AbstractParticleCollection) = isempty(get_particles(pc))

# to implement:
# Base.copy(pc::AbstractParticleCollection) = AbstractParticleCollection(copy(pc.particles))
# Base.empty(pc::AbstractParticleCollection) = AbstractParticleCollection(empty(pc.particles))
# Base.:(==)(pc1::AbstractParticleCollection, pc2::AbstractParticleCollection) = pc1.particles == pc2.particles

"""
    redshift(pc::AbstractParticleCollection) = 0

Returns the redshift at which the particles are located.

By default, 0 is returned. This method should be overridden for concrete types of `AbstractParticleCollection`
for which the redshift differs.
"""
redshift(pc::AbstractParticleCollection) = 0

"""
    CosmoParticles.show_properties(io::IO, mime, p::AbstractParticleCollection)

Prints the properties of the particles in the collection.

Should be used internally when overriding `Base.show` for concrete implementations of
`AbstractParticleCollection`.

This is not exported.
"""
show_properties(io::IO, mime::AbstractString, p::AbstractParticleCollection) = show_properties(io, MIME(mime), p)

function show_properties(io::IO, mime::MIME"text/plain", pc::AbstractParticleCollection)
    if !isempty(pc.particles)
        k = keys(pc) |> collect |> sort
        for key in k
            show(io, mime, pc[key])
            key != k[end] && println(io)
        end
    end
end
