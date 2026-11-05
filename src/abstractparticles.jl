"""
    abstract type AbstractParticles end

Abstract supertype for storing particle data efficiently in arrays.

Particles generally have IDs, positions, and velocities as properties.
The properties are expected to be stored in an `AbstractDict{Symbol,Any}`, such as the
[`OrderedDict`](@extref OrderedCollections.OrderedDict) of [`Particles`](@ref), as vectors of length ``N`` or
matrices of size ``m×N``, where ``N`` is the number of particles. Setting a property with a different number of
particles throws a `DimensionMismatch`. A scalar value, e.g. a mass that is equal for all particles, is stored as a
[`Fill`](https://github.com/JuliaArrays/FillArrays.jl) vector of length ``N``, so the particles must already have an
array property.

The inbuilt functionality of `AbstractParticles` includes accessing and setting properties
via the following syntax:
```jldoctest
julia> p = Particles(:dm);

julia> p.id = [1, 2, 3];

julia> p[:id] === p.id
true
```

A vector of property names selects these properties: the new particles share the arrays with `p`. A property is
removed with `delete!`.
```jldoctest
julia> p = Particles(:dm, :id => [1, 2, 3], :mass => [1.0, 2.0, 3.0], :temp => [10.0, 20.0, 30.0]);

julia> p2 = p[[:id, :mass]]
dm: 3 Particles
 id mass

julia> p2.id === p.id
true

julia> delete!(p, :temp)
dm: 3 Particles
 id mass
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
[`NamedRowArray`s](@extref NamedRowArrays.NamedRowArray).
```jldoctest
julia> p = Particles(:dm);

julia> # access by number always works, and by x, y, z for 2- and 3-dimensional matrices
       p.pos = [1.0 2.0; 3.0 4.0; 5.0 6.0];

julia> p.pos_1 == p.pos_x == @view p.pos[1, :]
true

julia> p.pos_3 == p.pos_z == @view p.pos[3, :]
true

julia> p.pos_zx # two of x, y, z return both rows in that order
2×2 view(::Matrix{Float64}, [3, 1], :) with eltype Float64:
 5.0  6.0
 1.0  2.0
```

```jldoctest
julia> using NamedRowArrays

julia> p = Particles(:dm, :arr => NamedRowArray([1 2; 3 4; 5 6; 7 8], [:a, :b, :c, :d]));

julia> p.arr_1 == p.arr_a == @view p.arr[1, :]
true

julia> p.arr_d
2-element view(::Matrix{Int64}, 4, :) with eltype Int64:
 7
 8
```

# Methods
The methods [`Base.keys`](@extref), [`Base.values`](@extref), [`Base.haskey`](@extref),
[`Base.empty!`](@extref), [`Base.isempty`](@extref), [`Base.copy!`](@extref), and [`Base.delete!`](@extref) with a
property name are forwarded to the property dictionary.

Concrete types of `AbstractParticles` should have the following methods implemented
(also see the implementation of [`Particles`](@ref)):
- [`Base.copy`](@extref): returns new object containing a copy of the dictionary
- [`Base.empty`](@extref): returns an identical object, but with an empty dictionary
- [`Base.:(==)`](@extref): should call [`Base.isequal`](@extref) to prevent obtaining `missing` results for
  properties with missing values (relevant for [`AllParticles`](@ref))
- [`Base.hash`](@extref): equal objects must have equal hashes
- [`CosmoParticles.particle_name`](@ref): returns the name of the struct to be printed via
  [`Base.show`](@extref Base.show-Tuple{IO, Any, Any})
- [`Base.propertynames`](@extref): implement this if there are additional struct fields
- [`Base.show(io, mime, p)`](@extref Base.show-Tuple{IO, Any, Any})

# Examples
```jldoctest
julia> using FillArrays

julia> p = Particles(:dm, :id => [1, 2, 3], :mass => Fill(2.0, 3));

julia> keys(p)
KeySet for a OrderedCollections.OrderedDict{Symbol, Any} with 2 entries. Keys:
  :id
  :mass
```
"""
abstract type AbstractParticles end

"""
    CosmoParticles.get_props(p::AbstractParticles)

Return the property dictionary belonging to the particles.

This returns `p.props` by default if not overridden.

This is meant for implementations of particle types.

# Examples
```jldoctest
julia> p = Particles(:dm, :id => [1, 2, 3]);

julia> CosmoParticles.get_props(p)
OrderedCollections.OrderedDict{Symbol, Any} with 1 entry:
  :id => [1, 2, 3]
```
"""
get_props(p::AbstractParticles) = p.props

# An `OrderedDict` of the requested type is stored as it is, so that the object shares it. Other dictionaries are
# copied into a new `OrderedDict`, since an unordered dictionary cannot be converted into one.
_ordered_dict(::Type{K}, ::Type{V}, d::OrderedDict{K,V}) where {K,V} = d
_ordered_dict(::Type{K}, ::Type{V}, d::AbstractDict) where {K,V} = OrderedDict{K,V}(d)

# All properties are arrays whose last dimension has the length of the particle number. Scalars are replaced by `Fill`
# vectors of that length.
function _check_props!(props::OrderedDict)
    n = nothing
    nkey = nothing
    for (key, val) in props
        val isa AbstractArray || continue
        m = size(val, ndims(val))
        if isnothing(n)
            n, nkey = m, key
        else
            m == n || throw(DimensionMismatch("property :$key has $m particles, but property :$nkey has $n"))
        end
    end
    for key in keys(props)
        props[key] isa AbstractArray && continue
        isnothing(n) && throw(_scalar_error(key))
        props[key] = Fill(props[key], n)
    end
    return props
end

_scalar_error(sym) = ArgumentError(
    "the scalar value of :$sym cannot be stored, since the particles have no array property that sets the " *
    "particle number",
)

# The particle number from the properties other than `sym`, or `nothing` if there are none.
function _particle_number_except(p::AbstractParticles, sym::Symbol)
    for (key, val) in get_props(p)
        key === sym && continue
        return size(val, ndims(val))
    end
    return nothing
end

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
    haskey(props, sym) && return props[sym]
    row = _get_matrix_row(p, sym)
    isnothing(row) || return row
    return getindex(props, sym)
end

# The underscore access `p.pos_x` returns a view of a row of the matrix property `pos`: by its number, by its name
# for a `NamedRowMatrix`, or by `x`, `y`, `z` for a matrix with two or three rows. Two axes, such as `p.pos_zx`,
# return a view of both rows in that order. Returns `nothing` if `sym` names no matrix row.
function _get_matrix_row(p::AbstractParticles, sym::Symbol)
    symsplit = split(String(sym), '_')
    for i in 1:(length(symsplit) - 1)
        key = Symbol(join(symsplit[1:i], "_"))
        key in keys(p) || continue

        vals = p[key]
        vals isa AbstractMatrix || continue
        row = join(symsplit[(i + 1):end], "_")
        row_sym = Symbol(row)
        row_int = tryparse(Int, row)

        !isnothing(row_int) && 1 ≤ row_int ≤ size(vals, 1) && return @view vals[row_int, :]
        vals isa NamedRowMatrix && row_sym in rownames(vals) && return @view vals[row_sym]
        if 2 ≤ size(vals, 1) ≤ 3 && row_sym in (:x, :y, :z)
            row_int = findfirst(==(row_sym), (:x, :y, :z))
            row_int ≤ size(vals, 1) && return @view vals[row_int, :]
        end
        if row_sym in (:xy, :yx) && 2 ≤ size(vals, 1) ≤ 3 ||
           row_sym in (:xz, :zx, :yz, :zy) && size(vals, 1) == 3
            return @view vals[get_dims(row_sym), :]
        end
    end
    return nothing
end

function Base.setindex!(p::AbstractParticles, val, sym::Symbol)
    n = _particle_number_except(p, sym)
    if val isa AbstractArray
        m = size(val, ndims(val))
        isnothing(n) ||
            m == n ||
            throw(DimensionMismatch("property :$sym has $m particles, but the particles have $n"))
    else
        isnothing(n) && throw(_scalar_error(sym))
        val = Fill(val, n)
    end
    return _setprop!(p, sym, val)
end

Base.delete!(p::AbstractParticles, sym::Symbol) = (delete!(get_props(p), sym); p)

# The selected properties share their arrays with `p`.
function Base.getindex(p::AbstractParticles, props::AbstractVector{Symbol})
    pnew = empty(p)
    for key in props
        key::Symbol # the elements of an `AbstractVector{Symbol}` are not inferred as `Symbol`
        _setprop!(pnew, key, p[key])
    end
    return pnew
end

# Stores a property without checks, for operations that keep the particle numbers of all properties consistent.
_setprop!(p::AbstractParticles, sym::Symbol, val) = (setindex!(get_props(p), val, sym); p)
Base.keys(p::AbstractParticles) = keys(get_props(p))
Base.haskey(p::AbstractParticles, key) =
    key in keys(p) || (key isa Symbol && !isnothing(_get_matrix_row(p, key)))

Base.values(p::AbstractParticles) = values(get_props(p))
Base.propertynames(p::AbstractParticles) = keys(p) |> collect
Base.empty!(p::AbstractParticles) = (empty!(get_props(p)); p)
Base.isempty(p::AbstractParticles) = isempty(get_props(p))
Base.copy!(dst::AbstractParticles, src::AbstractParticles) = (copy!(get_props(dst), get_props(src)); dst)

Base.getindex(p::AbstractParticles, ind::AbstractVector) = applyind(p, ind)
Base.deleteat!(p::AbstractParticles, ind::AbstractVector) = removeind!(p, ind)
"""
    CosmoParticles.deleteat(p::AbstractParticles, ind::AbstractVector)

Return a copy of the particles without the particles at the indices `ind`, leaving `p` unchanged.

This is an alias of [`CosmoParticles.removeind`](@ref), and the copying counterpart of
[`Base.deleteat!`](@extref).

# Examples
```jldoctest
julia> using FillArrays

julia> p = Particles(:dm, :id => [1, 2, 3], :mass => Fill(2.0, 3));

julia> CosmoParticles.deleteat(p, [2]).id
2-element Vector{Int64}:
 1
 3

julia> p.id
3-element Vector{Int64}:
 1
 2
 3
```
"""
deleteat(p::AbstractParticles, ind::AbstractVector) = removeind(p, ind)

# By default, the result has the properties of all particles, in the order in which they first appear. Particles
# without a property get `missing` values for it.
function Base.vcat(
    p::AbstractParticles,
    ps::AbstractParticles...;
    affect=unique(Iterators.flatten(keys.((p, ps...)))),
)
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
        _setprop!(p, key, dims == 1 ? reduce(vcat, vals) : reduce(hcat, vals))
    end

    return p
end

function _get_prop_dims(prop::Symbol, ps::AbstractParticles...)
    for p in ps
        vals = get(get_props(p), prop, missing)
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
[`Base.show`](@extref Base.show-Tuple{IO, Any, Any}).

# Examples
```jldoctest
julia> CosmoParticles.particle_name(Particles(:dm))
"Particles"
```
"""
function particle_name end

"""
    CosmoParticles.particle_number(p::AbstractParticles)

Returns the number of particles by determining the length of one of the property arrays.

# Examples
```jldoctest
julia> p = Particles(:dm, :id => [1, 2, 3], :pos => zeros(3, 3));

julia> CosmoParticles.particle_number(p)
3
```
"""
function particle_number(p::AbstractParticles)
    isempty(p) && return 0
    val = first(values(p))
    return size(val, ndims(val))
end

"""
    CosmoParticles.show_properties(io::IO, mime, p::AbstractParticles)

Prints the number of particles, the name of the type, and the property names.

Should be used internally when overriding [`Base.show`](@extref Base.show-Tuple{IO, Any, Any}) for concrete
implementations of [`AbstractParticles`](@ref).

# Examples
```jldoctest
julia> using FillArrays

julia> p = Particles(:dm, :id => [1, 2, 3], :mass => Fill(2.0, 3));

julia> CosmoParticles.show_properties(stdout, "text/plain", p)
3 Particles
 id mass
```
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
