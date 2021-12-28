struct ParticleCollection{AP<:AbstractParticles} <: AbstractParticleCollection{AP}
    particles::Dict{Symbol,AP}
end

ParticleCollection(::Type{AP}=Particles) where {AP<:AbstractParticles} = ParticleCollection{AP}(Dict{Symbol,AP}())
function ParticleCollection(pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles}
    ParticleCollection{AP}(Dict{Symbol,AP}(pairs...))
end

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


struct RedshiftParticleCollection{AP<:AbstractParticles,T<:Real} <: AbstractParticleCollection{AP}
    z::T
    particles::Dict{Symbol,AP}
end

function RedshiftParticleCollection(::Type{AP}, z::T) where {AP<:AbstractParticles,T<:Real}
    return ParticleCollection{AP,T}(z, Dict{Symbol,AP}())
end
function RedshiftParticleCollection(z::T) where {T<:Real}
    return ParticleCollection(Particles, z)
end
function RedshiftParticleCollection(z, pairs::Pair{Symbol,AP}...) where {AP<:AbstractParticles,T<:Real}
    RedshiftParticleCollection{AP,T}(z, Dict{Symbol,T}(pairs...))
end

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

redshift(pc::RedshiftParticleCollection) = pc.z

function Base.show(io::IO, mime::MIME"text/plain", pc::RedshiftParticleCollection)
    printstyled(io, "RedshiftParticleCollection"; bold=true)
    println(io, " at z = $(round(redshift(pc); sigdigits=4))")
    show_properties(io, mime, pc)
end
