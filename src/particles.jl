"""
    struct Particles <: AbstractParticles
        type::Symbol
        props::Dict{Symbol,Any}
    end

Particles of a certain type (typically something like `:dm` or `:gas` in the cosmological
context) with their properties.
"""
struct Particles <: AbstractParticles
    type::Symbol
    props::Dict{Symbol,Any}
end

"""
    Particles(type[, pairs::Pair...])

Create a [`Particles`](@ref) object with the given `type` and pairs of `Symbol` to the property
values. These are passed to the underlying `Dict`.
If no pairs are passed to the method, an empty `Dict` is created.
"""
Particles(type) = Particles(type, Dict{Symbol,Any}())
Particles(type, pairs::Pair...) = Particles(type, Dict{Symbol,Any}(pairs...))


Base.copy(p::Particles) = Particles(p.type, copy(p.props))
Base.:(==)(p1::Particles, p2::Particles) = p1.type == p2.type && p1.props == p2.props
particle_name(::Particles) = "Particles"

function Base.propertynames(p::Particles; private=false)
    if private
        return [fieldnames(Particles) |> collect; keys(p) |> collect]
    else
        return [:type; keys(p) |> collect]
    end
end

function Base.show(io::IO, mime::MIME"text/plain", p::Particles)
    printstyled(io, String(p.type); bold=true)
    print(io, ": ")
    show_properties(io, mime, p)
end
