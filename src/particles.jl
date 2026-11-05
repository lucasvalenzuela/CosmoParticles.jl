"""
    struct Particles <: AbstractParticles
        type::Symbol
        props::OrderedDict{Symbol,Any}
    end

Particles of a certain type (typically something like `:dm` or `:gas` in the cosmological
context) with their properties.

The following naming convention is to be used for specific properties:
- `:id`: vector of IDs
- `:pos`: ``2×N`` or ``3×N`` matrix with positions
- `:vel`: ``2×N`` or ``3×N`` matrix with velocities
- `:mass`: vector of masses
- `:temp`: vector of temperatures

# Examples
```jldoctest
julia> p = Particles(:gas, :id => [1, 2, 3], :temp => [1e4, 2e4, 3e4])
gas: 3 Particles
 id temp

julia> p.type
:gas

julia> p[:temp] # equivalent to p.temp
3-element Vector{Float64}:
 10000.0
 20000.0
 30000.0
```
"""
struct Particles <: AbstractParticles
    type::Symbol
    props::OrderedDict{Symbol,Any}

    Particles(type::Symbol, props::AbstractDict) = new(type, _check_props!(_ordered_dict(Symbol, Any, props)))
end

"""
    Particles(type::Symbol[, pairs::Pair...])

Create a [`Particles`](@ref) object with the given `type` and pairs of `Symbol` to the property
values. These are passed to the underlying [`OrderedDict`](@extref OrderedCollections.OrderedDict), which keeps the
order of the properties. If no pairs are passed to the method, an empty `OrderedDict` is created.

The properties can also be passed as an `AbstractDict`, as in `Particles(type, props)`. An `OrderedDict{Symbol,Any}`
is stored as it is, and other dictionaries are copied into a new `OrderedDict`.

# Examples
```jldoctest
julia> using FillArrays

julia> p = Particles(:dm, :id => [1, 2, 3], :mass => Fill(2.0, 3));

julia> p.mass
3-element Fill{Float64}, with entries equal to 2.0

julia> Particles(:dm, Dict(:id => [1, 2, 3])).id
3-element Vector{Int64}:
 1
 2
 3
```
"""
Particles(type::Symbol) = Particles(type, OrderedDict{Symbol,Any}())
Particles(type::Symbol, pairs::Pair...) = Particles(type, OrderedDict{Symbol,Any}(pairs...))


Base.copy(p::Particles) = Particles(p.type, copy(p.props))
Base.empty(p::Particles) = Particles(p.type, empty(p.props))
Base.:(==)(p1::Particles, p2::Particles) = p1.type == p2.type && isequal(p1.props, p2.props)
Base.hash(p::Particles, h::UInt) = hash(p.props, hash(p.type, hash(:Particles, h)))
particle_name(::Particles) = "Particles"

function Base.propertynames(p::Particles, private::Bool=false)
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
