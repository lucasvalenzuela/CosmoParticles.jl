"""
    abstract type AbstractParticles end

Abstract supertype for storing particle data efficiently in arrays.

Particles generally have IDs, positions, and velocities as properties.
The properties are expected to be stored in a `Dict{Symbol,Any}` as vectors of length ``N``,
matrices of size ``m×N``, or as scalars (when the property is equal for all particles).

The inbuilt functionality of `AbstractParticles` includes accessing and setting properties
via the following syntax:
```julia
p.id = [1, 2, 3]
p[:id] === p.id
```

If the struct has additional fields, these can also be accessed by `p.field`, but not by
`p[:field]`. The latter syntax can only be used to access the particle properties.

# Property keys
- `:id`: vector of IDs
- `:pos`: ``2×N`` or ``3×N`` matrix with positions
- `:vel`: ``2×N`` or ``3×N`` matrix with velocities

Any arrays may be of the type `AbstractArray`, provided the arrays are 1-indexed.

# Methods
The methods `Base.keys`, `Base.values`, `Base.haskey`, `Base.empty`, `Base.empty!`, `Base.isempty`,
and `Base.copy!` are forwarded to the property `Dict`.

Concrete types of `AbstractParticles` should have the following methods implemented
(also see the implementation of [`Particles`](@ref)):
- `Base.copy`: returns new object containing a copy of the `Dict`
- `Base.empty`: returns an identical object, but with an empty `Dict`
- `Base.:(==)`: should call `Base.isequal` to prevent obtaining `missing` results for properties with missing values (relevant for [`AllParticles`](@ref))
- [`CosmoParticles.particle_name`](@ref): returns the name of the struct to be printed via `Base.show`
- `Base.propertynames`: implement this if there are additional struct fields
- `Base.show(io, mime, p)`
"""
abstract type AbstractParticles end

"""
    CosmoParticles.get_props(p::AbstractParticles)

Return the property `Dict` belonging to the particles.

This returns `p.props` by default if not overridden.

This is not exported and should not be used outside of the defining files for particle types.
"""
get_props(p::AbstractParticles) = p.props

function Base.getproperty(p::AP, sym::Symbol) where {AP<:AbstractParticles}
    if sym in fieldnames(AP)
        return getfield(p, sym)
    else
        return getindex(p, sym)
    end
end

function Base.setproperty!(p::AP, sym::Symbol, val) where {AP<:AbstractParticles}
    if sym in fieldnames(AP)
        setfield!(p, sym, val)
    else
        setindex!(p, val, sym)
    end
end

Base.getindex(p::AbstractParticles, sym::Symbol) = getindex(get_props(p), sym)
Base.setindex!(p::AbstractParticles, val, sym::Symbol) = (setindex!(get_props(p), val, sym); p)
Base.keys(p::AbstractParticles) = keys(get_props(p))
Base.haskey(p::AbstractParticles, key) = key in keys(p)
Base.values(p::AbstractParticles) = values(get_props(p))
Base.propertynames(p::AbstractParticles) = keys(p) |> collect
Base.empty!(p::AbstractParticles) = (empty!(get_props(p)); p)
Base.isempty(p::AbstractParticles) = isempty(get_props(p))
Base.copy!(dst::AbstractParticles, src::AbstractParticles) = (copy!(get_props(dst), get_props(src)); dst)

Base.getindex(p::AbstractParticles, ind::AbstractVector) = applyind(p, ind)
Base.deleteat!(p::AbstractParticles, ind::AbstractVector) = removeind!(p, ind)
deleteat(p::AbstractParticles, ind::AbstractVector) = removeind(p, ind)

# to implement:
# Base.copy(p::AbstractParticles) = AbstractParticles(copy(p.props))
# Base.empty(p::AbstractParticles) = AbstractParticles(empty(p.props))
# Base.:(==)(p1::AbstractParticles, p2::AbstractParticles) = isequal(p1.props, p2.props)

function Base.vcat(p::AbstractParticles, ps::AbstractParticles...; affect=keys(p))
    pout = empty(p)

    for key in affect
        dims = _get_prop_dims(key, p, ps...)
        vals = _get_prop_scaled.(key, dims, [p, ps...])
        pout[key] = dims == 1 ? reduce(vcat, vals) : reduce(hcat, vals)
    end

    return pout
end

function Base.append!(p::AbstractParticles, ps::AbstractParticles...)
    n = particle_number.([p, ps...])
    for key in keys(p)
        dims = _get_prop_dims(key, p, ps...)
        vals = _get_prop_scaled.(key, dims, [p, ps...], n)
        p[key] = dims == 1 ? reduce(vcat, vals) : reduce(hcat, vals)
    end

    return p
end

function _get_prop_dims(prop::Symbol, ps::AbstractParticles...)
    for p in ps
        vals = p[prop]
        if vals isa AbstractVector
            return 1
        elseif vals isa AbstractMatrix
            return size(vals, 1)
        end
    end

    return 1
end

function _get_prop_scaled(prop::Symbol, dims::Integer, p::AbstractParticles, n=particle_number(p))
    vals = get(get_props(p), prop, missing)

    if vals isa AbstractVecOrMat
        return vals
    else
        if dims == 1
            return fill(vals, n)
        else
            return fill(vals, dims, n)
        end
    end
end

"""
    CosmoParticles.particle_name(p::AbstractParticles)

Returns the name of the particle type, which is used when printing an object of this type via
`Base.show`.

This is not exported.
"""
function particle_name(::AbstractParticles) end

"""
    CosmoParticles.particle_number(p::AbstractParticles)

Returns the number of particles by determining the length of one of the property arrays.

This is not exported.
"""
function particle_number(p::AbstractParticles)
    for val in values(p)
        if val isa AbstractArray
            return size(val, ndims(val))
        end
    end
    return 0
end

"""
    CosmoParticles.show_properties(io::IO, mime, p::AbstractParticles)

Prints the number of particles, the name of the type, and the property names.

Should be used internally when overriding `Base.show` for concrete implementations of
`AbstractParticles`.

This is not exported.
"""
show_properties(io::IO, mime::AbstractString, p::AbstractParticles) = show_properties(io, MIME(mime), p)

function show_properties(io::IO, ::MIME"text/plain", p::AbstractParticles)
    n = particle_number(p)
    name = particle_name(p)
    println(io, "$n $name")
    print(io, " ")
    props = join(keys(p) |> collect |> sort, " ")
    isempty(props) || print(io, props)
end
