"""
    AllParticles{APC} <: AbstractParticles where {APC<:AbstractParticleCollection}

Combines the particles of different types from a particle collection lazily.

The property arrays of the collection's particles are concatenated lazily with
[LazyArrays.jl](https://github.com/JuliaArrays/LazyArrays.jl).
Scalar properties (e.g. mass for equal-mass particles) are repeated according to the particle number of the given
particle type with [FillArrays.jl](https://github.com/JuliaArrays/FillArrays.jl).
Properties that only exist for one of the particle types are repeated as `missing` values in the same way as
scalar properties.

# Usage
The individual properties of `AllParticles` can be converted as regular arrays with `Array(ap.pos)`
and all or selected particle properties can be converted to a `Particles` object with `Particles(ap)`,
where all lazy arrays are materialized as `Array`s.
`AllParticles` cannot be sorted, filtered, or copied.
In-place modifications are possible, such as [`rotate!`](@ref), [`translate!`](@ref), or [`to_comoving`](@ref).

!!! warning "In-place modifications"
    Modifying an `AllParticles` object in-place results in the modification of the underlying particle collection!
    To prevent bugs or unexpected behavior, it is advised to modify the particle collection instead and then
    obtain the object through `pc.all` or `AllParticles(pc)` (these are equivalent).

The lazy concatenation occurs only when the property is needed and is only performed once.
Note that this means that calling `pc.all.prop` multiple times results in the given property being lazily
concatenated every time. As long as the property array pointers of the particle collection are not modified
(e.g. by `sort!` or `filter!`), it is possible to reuse an object `ap = pc.all` without the need of regenerating
the concatenated arrays.

!!! warning "Modifying the underlying particle collection"
    After modifying the property array pointers of a particle collection (e.g. by `sort!` or `filter!`), any
    already created `AllParticles` objects from that collection may have lazy arrays pointing to the old,
    non-modified arrays. When this happens, create a new object through `pc.all` or `AllParticles(pc)`.

# Example
```julia
julia> pc::ParticleCollection
ParticleCollection
dm: 100 Particles
 id mass pos
gas: 50 Particles
 id mass pos temp

julia> ap = AllParticles(pc) # equivalent to ap = pc.all
all: 150 Particles
 id mass pos temp

julia> p = Particles(ap) # p is detached from pc and therefore modifiable
all: 150 Particles
 id mass pos temp

julia> sort!(p, :id)
```
"""
struct AllParticles{APC} <: AbstractParticles where {APC<:AbstractParticleCollection}
    particle_collection::APC
    props::Dict{Symbol,Any}

    AllParticles(pc::APC) where {APC} = new{APC}(pc, Dict{Symbol,Any}())
end


"""
    particle_collection(p::AllParticles)

Returns the underlying particle collection.

This is not exported.
"""
particle_collection(p::AllParticles) = p.particle_collection


function Base.getindex(ap::AllParticles, sym::Symbol)
    # get property index if property already saved or if key does not exist at all to raise Dict error
    if haskey(ap.props, sym) || !haskey(ap, sym)
        return ap.props[sym]
    end

    pc = ap.particle_collection

    # get array dimensions
    dims = 0 # 0 for vectors; for arrays: size(arr, 1)
    for key in keys(pc)
        p = pc[key]
        if haskey(p, sym)
            p[sym] isa AbstractVector && break
            if p[sym] isa AbstractMatrix
                dims = size(p[sym], 1)
                break
            end
        end
    end

    # create array of particle property arrays
    arrs = []

    # iterate through particles
    for p in values(pc)
        # fill array with missing if property does not exist
        if haskey(p, sym)
            # fill array with scalar value if property is not an array
            if p[sym] isa AbstractArray
                arr = p[sym]
            else
                if dims == 0
                    arr = Fill(p[sym], particle_number(p))
                else
                    arr = Fill(p[sym], dims, particle_number(p))
                end
            end
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

function Base.keys(ap::AllParticles)
    keys_arr = [keys(p) for p in values(ap.particle_collection)]
    length(keys_arr) == 1 && return keys_arr[1]
    return union(keys_arr...) |> collect
end

Base.haskey(p::AllParticles, key) = haskey(p.props, key) || key in keys(p)

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

function Base.copy!(::AllParticles, ::AllParticles, props=())
    error(
        "Copying particles to an AllParticles object is not allowed. Try converting it via Particles(allparticles) instead.",
    )
end

function Base.copy!(::AllParticles, ::AbstractParticles, props=())
    error(
        "Copying particles to an AllParticles object is not allowed. Try converting it via Particles(allparticles) instead.",
    )
end

"""
    Base.copy!(dst::AbstractParticles, src::AllParticles[, props])

Copies the materialized particle properties to another particle object.

To only copy certain properties, `props` can be passed as a tuple of `Symbol`s.
The `dst` cannot be of type `AllParticles`.
"""
function Base.copy!(dst::AbstractParticles, src::AllParticles, props=())
    empty!(dst)
    for key in keys(src)
        if isempty(props) || key in props
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
    for p in values(p.particle_collection)
        isempty(p) || return false
    end
    return true
end

Base.:(==)(p1::AllParticles, p2::AllParticles) = p1.particle_collection == p2.particle_collection

particle_name(::AllParticles) = "Particles"

function particle_number(ap::AllParticles)
    n = 0
    for p in values(ap.particle_collection)
        n += particle_number(p)
    end
    return n
end

function Base.propertynames(p::AllParticles; private=false)
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
    Particles(p::AllParticles, props=())

Materializes the particle properties in a new particle object of type `:all`.

To only copy certain properties, `props` can be passed as a tuple of `Symbol`s.
"""
Particles(p::AllParticles, props=()) = copy!(Particles(:all), p, props)

function applyind!(::AllParticles, ::AbstractVector)
    error(
        "Indexing into an AllParticles object is not allowed. Try converting it via Particles(allparticles) first.",
    )
end

applyind(p::AllParticles, ind::AbstractVector) = applyind!(p, ind)
