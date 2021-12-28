abstract type AbstractParticleCollection{AP} end

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
