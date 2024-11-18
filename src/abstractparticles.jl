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

For matrix properties, a simplified syntax with underscores can be used to access the matrix rows. These are
returned as views into the array, so `collect` them if necessary. This can also be used with the row names of
[`NamedRowArray`s](https://lucasvalenzuela.github.io/NamedRowArrays.jl/dev/).
```julia
# access by number always works, and by x, y, z for 2- and 3-dimensional matrices
p.pos = rand(3, 100)
p.pos_1 == p.pos_x == @view p.pos[1, :]
p.pos_2 == p.pos_y == @view p.pos[2, :]
p.pos_3 == p.pos_z == @view p.pos[3, :]

p.arr = NamedRowArray(rand(4, 100), [:a, :b, :c, :d])
p.arr_1 == p.arr_a == @view p.arr[1, :] == @view p.arr[:a]
```

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

function Base.getindex(p::AbstractParticles, sym::Symbol)
    props = get_props(p)
    if !haskey(props, sym)
        symsplit = split(String(sym), '_')
        for i in 1:(length(symsplit) - 1)
            symtest = join(symsplit[1:i], "_") |> Symbol

            # key exists
            if haskey(props, symtest)
                vals = getindex(props, symtest)
                key = join(symsplit[i+1:end], "_")
                key_sym = Symbol(key)
                key_int = tryparse(Int, key)

                if vals isa AbstractMatrix && !isnothing(key_int)
                    1 ≤ key_int ≤ size(vals, 1) && return @view vals[key_int, :]
                end

                if vals isa NamedRowMatrix
                    key_sym in names(vals) && return @view vals[key_sym]
                end

                if vals isa AbstractMatrix && 2 ≤ size(vals, 1) ≤ 3 && key_sym in (:x, :y, :z)
                    key_int = findfirst(==(key_sym), (:x, :y, :z))
                    1 ≤ key_int ≤ size(vals, 1) && return @view vals[key_int, :]
                end
            end
        end
    end

    return getindex(props, sym)
end

Base.setindex!(p::AbstractParticles, val, sym::Symbol) = (setindex!(get_props(p), val, sym); p)
Base.keys(p::AbstractParticles) = keys(get_props(p))

function Base.haskey(p::AbstractParticles, key)
    if key in keys(p)
        return true
    end

    key isa Symbol || return false

    props = get_props(p)
    symsplit = split(String(key), '_')
    for i in 1:(length(symsplit) - 1)
        symtest = join(symsplit[1:i], "_") |> Symbol

        # key exists
        if haskey(props, symtest)
            vals = getindex(props, symtest)
            key = join(symsplit[i+1:end], "_")
            key_sym = Symbol(key)
            key_int = tryparse(Int, key)

            if vals isa AbstractMatrix && !isnothing(key_int)
                1 ≤ key_int ≤ size(vals, 1) && return true
            end

            if vals isa NamedRowMatrix
                key_sym in names(vals) && return true
            end

            if vals isa AbstractMatrix && 2 ≤ size(vals, 1) ≤ 3 && key_sym in (:x, :y, :z)
                key_int = findfirst(==(key_sym), (:x, :y, :z))
                1 ≤ key_int ≤ size(vals, 1) && return true
            end
        end
    end

    return false
end

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
