"""
    AllParticles{APC} <: AbstractParticles where {APC<:AbstractParticleCollection}

Combines the particles of different types from a particle collection lazily.

The property arrays of the collection's particles are concatenated lazily with
[LazyArrays.jl](https://github.com/JuliaArrays/LazyArrays.jl).
Properties that do not exist for a particle type are filled with `missing` values for its particles with
[FillArrays.jl](https://github.com/JuliaArrays/FillArrays.jl).

# Usage
The individual properties of `AllParticles` can be converted as regular arrays with `Array(ap.pos)`
and all or selected particle properties can be converted to a [`Particles`](@ref) object with `Particles(ap)`,
where all lazy arrays are materialized as `Array`s.
`AllParticles` cannot be sorted, filtered, or copied, and the in-place transformations [`rotate!`](@ref),
[`translate!`](@ref), [`translate_periodic!`](@ref), [`translate_periodic_to_center!`](@ref),
[`to_comoving!`](@ref), and [`to_physical!`](@ref) throw an error for it. They should be applied to the particle
collection instead.

# Creation
```jldoctest
julia> pc = ParticleCollection(
           Particles(:dm, :id => [1, 2]), Particles(:stars, :id => [3]), Particles(:gas, :id => [4]),
       );

julia> ap = AllParticles(pc);

julia> ap == pc.all
true

julia> ap = AllParticles(pc, [:dm, :stars]); # just include certain particle types

julia> Array(ap.id)
3-element Vector{Int64}:
 1
 2
 3
```

!!! warning "Writing into the lazy arrays"
    The lazy arrays refer to the arrays of the particle collection, so writing into them, e.g. with `ap.pos .= 0`,
    modifies the particle collection. To prevent bugs or unexpected behavior, modify the particle collection instead
    and then obtain the object through `pc.all` or `AllParticles(pc)` (these are equivalent).

The lazy concatenation occurs only when the property is needed and is only performed once.
Note that this means that calling `pc.all.prop` multiple times results in the given property being lazily
concatenated every time. As long as the property array pointers of the particle collection are not modified
(e.g. by [`sort!`](@ref Base.sort!) or [`filter!`](@ref Base.filter!)), it is possible to reuse an object `ap = pc.all`
without the need of regenerating the concatenated arrays.

!!! warning "Modifying the underlying particle collection"
    After modifying the property array pointers of a particle collection (e.g. by `sort!` or `filter!`), any
    already created `AllParticles` objects from that collection may have lazy arrays pointing to the old,
    non-modified arrays. When this happens, create a new object through `pc.all` or `AllParticles(pc)`.

# Examples
```jldoctest
julia> using FillArrays

julia> pc = ParticleCollection(
           Particles(:dm, :id => [3, 1], :mass => Fill(2.0, 2)),
           Particles(:gas, :id => [2], :mass => Fill(1.0, 1), :temp => [1e4]),
       );

julia> ap = AllParticles(pc) # equivalent to ap = pc.all
all: 3 Particles
 id mass temp

julia> p = Particles(ap) # p is detached from pc and therefore modifiable
all: 3 Particles
 id mass temp

julia> sort!(p, :id);

julia> p.id
3-element Vector{Int64}:
 1
 2
 3

julia> Particles(ap, [:id, :mass]) # only get certain properties
all: 3 Particles
 id mass
```
"""
struct AllParticles{APC<:AbstractParticleCollection} <: AbstractParticles
    particle_collection::APC
    ptypes::Vector{Symbol}
    props::OrderedDict{Symbol,Any}

    function AllParticles(pc::APC) where {APC}
        ptypes = collect(keys(pc))
        return new{APC}(pc, ptypes, OrderedDict{Symbol,Any}())
    end

    function AllParticles(pc::APC, ptypes) where {APC}
        if !(ptypes isa Vector)
            ptypes = collect(ptypes)
        end

        for key in ptypes
            haskey(pc, key) || throw(ArgumentError("the particle collection does not have the particle type :$key"))
        end

        return new{APC}(pc, ptypes, OrderedDict{Symbol,Any}())
    end
end


"""
    CosmoParticles.particle_collection(p::AllParticles)

Returns the underlying particle collection.

# Examples
```jldoctest
julia> pc = ParticleCollection(Particles(:dm, :id => [1, 2]));

julia> CosmoParticles.particle_collection(pc.all) === pc
true
```
"""
particle_collection(p::AllParticles) = p.particle_collection


# The concatenated array requires vectors for all particle types, or matrices with the same number of rows.
function _check_alldims(arr::AbstractArray, dims, sym, ptype, dims_ptype)
    arrdims = arr isa AbstractVector ? 0 : size(arr, 1)
    arrdims == dims && return nothing
    desc(d) = d == 0 ? "a vector" : "a matrix with $d rows"
    throw(
        DimensionMismatch(
            "property :$sym is $(desc(arrdims)) for the :$ptype particles, but $(desc(dims)) for the :$dims_ptype particles",
        ),
    )
end

function Base.getindex(ap::AllParticles, sym::Symbol)
    haskey(ap.props, sym) && return ap.props[sym]
    if !(sym in keys(ap))
        row = _get_matrix_row(ap, sym)
        isnothing(row) || return row
        return ap.props[sym] # throws a KeyError
    end

    pc = ap.particle_collection

    # get array dimensions
    dims = 0 # 0 for vectors; for arrays: size(arr, 1)
    dims_ptype = nothing # particle type that determines `dims`
    for key in ap.ptypes
        p = pc[key]
        if haskey(p, sym)
            if p[sym] isa AbstractVector
                dims_ptype = key
                break
            end
            if p[sym] isa AbstractMatrix
                dims = size(p[sym], 1)
                dims_ptype = key
                break
            end
        end
    end

    # create array of particle property arrays
    arrs = []

    # iterate through particles
    for key in ap.ptypes
        p = pc[key]
        # fill array with missing if property does not exist
        if haskey(p, sym)
            arr = p[sym]
            _check_alldims(arr, dims, sym, key, dims_ptype)
            push!(arrs, arr)
        else
            if dims == 0
                arr = Fill(missing, particle_number(p))
            else
                arr = Fill(missing, dims, particle_number(p))
            end
            push!(arrs, arr)
        end
    end

    cat_func = dims == 0 ? vcat : hcat
    lazy = ApplyArray(cat_func, arrs...)

    ap.props[sym] = lazy
    return lazy
end

function Base.setindex!(::AllParticles, _, ::Symbol)
    error(
        "Setting a property of an AllParticles object is not allowed. Try converting it via Particles(allparticles) first.",
    )
end

_setprop!(ap::AllParticles, sym::Symbol, val) = setindex!(ap, val, sym)

function Base.delete!(::AllParticles, ::Symbol)
    error(
        "Deleting a property of an AllParticles object is not allowed. Try converting it via Particles(allparticles) first.",
    )
end

function Base.getindex(::AllParticles, ::AbstractVector{Symbol})
    error(
        "Selecting properties of an AllParticles object is not supported. Try Particles(allparticles, props) instead.",
    )
end

function Base.keys(ap::AllParticles)
    keys_arr = [keys(ap.particle_collection[ptype]) for ptype in ap.ptypes]
    length(keys_arr) == 0 && return Symbol[]
    length(keys_arr) == 1 && return keys_arr[1]
    return unique(Iterators.flatten(keys_arr))
end

Base.values(p::AllParticles) = Iterators.map(key -> getindex(p, key), keys(p))

function Base.copy(::AllParticles)
    error(
        "Copying an AllParticles object is not allowed. Try converting it via Particles(allparticles) instead.",
    )
    ## NOTE: potential way of making copy work
    # pnew = Particles(p.particle_collection)
    # copy!(pnew.props, p.props)
    # return pnew
end

function Base.copy!(::AllParticles, ::AllParticles, props=nothing)
    error(
        "Copying particles to an AllParticles object is not allowed. Try converting it via Particles(allparticles) instead.",
    )
end

function Base.copy!(::AllParticles, ::AbstractParticles, props=nothing)
    error(
        "Copying particles to an AllParticles object is not allowed. Try converting it via Particles(allparticles) instead.",
    )
end

"""
    Base.copy!(dst::AbstractParticles, src::AllParticles[, props])

Copies the materialized particle properties to another particle object.

To only copy certain properties, `props` can be passed as a vector or tuple of `Symbol`s.
The `dst` cannot be of type [`AllParticles`](@ref).

# Examples
```jldoctest
julia> using FillArrays

julia> dm = Particles(:dm, :id => [1, 2], :mass => Fill(2.0, 2));

julia> gas = Particles(:gas, :id => [3], :mass => Fill(1.0, 1));

julia> pc = ParticleCollection(dm, gas);

julia> p = Particles(:all);

julia> copy!(p, pc.all, [:id]);

julia> p.id
3-element Vector{Int64}:
 1
 2
 3
```
"""
function Base.copy!(dst::AbstractParticles, src::AllParticles, props=nothing)
    empty!(dst)
    for key in keys(src)
        if isnothing(props) || key in props
            val = src[key]
            dst[key] = Array(val)
        end
    end
    return dst
end

function Base.empty(::AllParticles)
    error(
        "Creating an empty AllParticles object is not supported. Try emptying the underlying particle collection instead.",
    )
end

function Base.empty!(::AllParticles)
    error(
        "Emptying an AllParticles object is not allowed. Try emptying the underlying particle collection instead.",
    )
end

function Base.isempty(p::AllParticles)
    pc = p.particle_collection
    for ptype in p.ptypes
        isempty(pc[ptype]) || return false
    end
    return true
end

Base.:(==)(p1::AllParticles, p2::AllParticles) =
    p1.particle_collection == p2.particle_collection && p1.ptypes == p2.ptypes
Base.hash(p::AllParticles, h::UInt) = hash(p.ptypes, hash(p.particle_collection, hash(:AllParticles, h)))

particle_name(::AllParticles) = "Particles"

function particle_number(ap::AllParticles)
    n = 0
    pc = ap.particle_collection
    for ptype in ap.ptypes
        n += particle_number(pc[ptype])
    end
    return n
end

function Base.propertynames(p::AllParticles, private::Bool=false)
    if private
        return [fieldnames(AllParticles) |> collect; keys(p) |> collect]
    else
        return keys(p) |> collect
    end
end

function Base.show(io::IO, mime::MIME"text/plain", p::AllParticles)
    printstyled(io, "all"; bold=true)
    print(io, ": ")
    show_properties(io, mime, p)
end


"""
    Particles(p::AllParticles, props=nothing)

Materializes the particle properties in a new particle object of type `:all`.

To only copy certain properties, `props` can be passed as a vector or tuple of `Symbol`s.

# Examples
```jldoctest
julia> pc = ParticleCollection(Particles(:dm, :id => [1, 2]), Particles(:gas, :id => [3]));

julia> p = Particles(pc.all);

julia> p.id
3-element Vector{Int64}:
 1
 2
 3
```
"""
Particles(p::AllParticles, props=nothing) = copy!(Particles(:all), p, props)

function applyind!(::AllParticles, ::AbstractVector)
    error(
        "Indexing into an AllParticles object is not allowed. Try converting it via Particles(allparticles) first.",
    )
end

applyind(p::AllParticles, ind::AbstractVector) = applyind!(p, ind)
