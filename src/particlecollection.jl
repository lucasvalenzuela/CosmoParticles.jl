"""
    struct ParticleCollection{AP<:AbstractParticles} <: AbstractParticleCollection{AP}
        particles::OrderedDict{Symbol,AP}
    end

Basic particle collection at redshift 0.

`copy(pc)` copies the dictionary of the collection, but not the particles: the copy contains the same
[`Particles`](@ref) objects as `pc`.

# Examples
```jldoctest
julia> pc = ParticleCollection(Particles(:dm, :id => [1, 2]), Particles(:gas, :id => [3]))
ParticleCollection
dm: 2 Particles
 id
gas: 1 Particles
 id

julia> pc[:gas] # equivalent to pc.gas
gas: 1 Particles
 id

julia> pc.gas.id
1-element Vector{Int64}:
 3
```
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

An `OrderedDict{Symbol,AP}` is stored as it is, and other dictionaries are copied into a new
[`OrderedDict`](@extref OrderedCollections.OrderedDict).

# Examples
```jldoctest
julia> d = Dict(:dm => Particles(:dm, :id => [1, 2]));

julia> keys(ParticleCollection(d))
KeySet for a OrderedCollections.OrderedDict{Symbol, Particles} with 1 entry. Keys:
  :dm
```
"""
function ParticleCollection(particles::AbstractDict{Symbol,AP}) where {AP<:AbstractParticles}
    return ParticleCollection{AP}(particles)
end

# Each particle type can appear only once in a collection.
function _particles_dict(::Type{AP}, pairs) where {AP<:AbstractParticles}
    particles = OrderedDict{Symbol,AP}()
    for (ptype, p) in pairs
        haskey(particles, ptype) && throw(ArgumentError("the particle type :$ptype is given more than once"))
        particles[ptype] = p
    end
    return particles
end

"""
    ParticleCollection(::Type{AP}=Particles) where {AP<:AbstractParticles}

Create an empty particle collection of the passed type.

By default, the particle collection is created for the [`Particles`](@ref) type.

# Examples
```jldoctest
julia> pc = ParticleCollection();

julia> isempty(pc)
true
```
"""
ParticleCollection(::Type{AP}=Particles) where {AP<:AbstractParticles} =
    ParticleCollection{AP}(OrderedDict{Symbol,AP}())

"""
    ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}

Create a particle collection from pairs of Symbols and particles.

Each particle type can only be given once.

# Examples
```jldoctest
julia> pc = ParticleCollection(:dm => Particles(:dm), :gas => Particles(:gas));

julia> collect(keys(pc))
2-element Vector{Symbol}:
 :dm
 :gas
```
"""
function ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}
    ParticleCollection{AP}(_particles_dict(AP, pairs))
end

"""
    ParticleCollection(p::Particles...)

Create a particle collection from multiple [`Particles`](@ref).

Each particle type can only be given once.

The keys are taken from the [`Particles`](@ref) type.

# Examples
```jldoctest
julia> pc = ParticleCollection(Particles(:dm), Particles(:gas));

julia> collect(keys(pc))
2-element Vector{Symbol}:
 :dm
 :gas
```
"""
function ParticleCollection(p::Particles, ps::Particles...)
    return ParticleCollection(
        _particles_dict(Particles, (particle.type => particle for particle in (p, ps...))),
    )
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

`copy(pc)` copies the dictionary of the collection, but not the particles: the copy contains the same
[`Particles`](@ref) objects as `pc`.

# Examples
```jldoctest
julia> pc = RedshiftParticleCollection(0.5, Particles(:dm, :id => [1, 2]))
ParticleCollection at z = 0.5
dm: 2 Particles
 id

julia> redshift(pc)
0.5
```
"""
struct RedshiftParticleCollection{AP<:AbstractParticles,T<:Real} <: AbstractParticleCollection{AP}
    z::T
    particles::OrderedDict{Symbol,AP}

    function RedshiftParticleCollection{AP,T}(
        z,
        particles::AbstractDict,
    ) where {AP<:AbstractParticles,T<:Real}
        return new{AP,T}(z, _ordered_dict(Symbol, AP, particles))
    end
end

"""
    RedshiftParticleCollection(z::T, particles::AbstractDict{Symbol,AP}) where {AP<:AbstractParticles,T<:Real}

Create a particle collection at redshift `z` from a dictionary of particles.

An `OrderedDict{Symbol,AP}` is stored as it is, and other dictionaries are copied into a new
[`OrderedDict`](@extref OrderedCollections.OrderedDict).

# Examples
```jldoctest
julia> d = Dict(:dm => Particles(:dm, :id => [1, 2]));

julia> redshift(RedshiftParticleCollection(0.5, d))
0.5
```
"""
function RedshiftParticleCollection(
    z::T,
    particles::AbstractDict{Symbol,AP},
) where {AP<:AbstractParticles,T<:Real}
    return RedshiftParticleCollection{AP,T}(z, particles)
end

"""
    RedshiftParticleCollection([::Type{AP}=Particles,] z::T) where {AP<:AbstractParticles,T<:Real}

Create an empty particle collection of the passed type at redshift `z`.

By default, the particle collection is created for the [`Particles`](@ref) type.

# Examples
```jldoctest
julia> pc = RedshiftParticleCollection(0.5);

julia> isempty(pc), redshift(pc)
(true, 0.5)
```
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

Each particle type can only be given once.

# Examples
```jldoctest
julia> pc = RedshiftParticleCollection(0.5, :dm => Particles(:dm), :gas => Particles(:gas));

julia> redshift(pc), collect(keys(pc))
(0.5, [:dm, :gas])
```
"""
function RedshiftParticleCollection(z::T, pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles,T<:Real}
    RedshiftParticleCollection{AP,T}(z, _particles_dict(AP, pairs))
end

"""
    RedshiftParticleCollection(z::Real, p::Particles...)

Create a particle collection at redshift `z` from multiple [`Particles`](@ref).

Each particle type can only be given once.

The keys are taken from the [`Particles`](@ref) type.

# Examples
```jldoctest
julia> pc = RedshiftParticleCollection(0.5, Particles(:dm), Particles(:gas));

julia> redshift(pc), collect(keys(pc))
(0.5, [:dm, :gas])
```
"""
function RedshiftParticleCollection(z::Real, p::Particles, ps::Particles...)
    particles = _particles_dict(Particles, (particle.type => particle for particle in (p, ps...)))
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

# A `ParticleCollection` equals a `RedshiftParticleCollection` at redshift 0, so the redshift enters the hash only
# when it is not zero.
Base.hash(pc::ParticleCollection, h::UInt) = hash(pc.particles, hash(:ParticleCollection, h))
function Base.hash(pc::RedshiftParticleCollection, h::UInt)
    h = hash(pc.particles, hash(:ParticleCollection, h))
    return iszero(pc.z) ? h : hash(pc.z, h)
end
