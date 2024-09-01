"""
    Base.sort!(p::AbstractParticles, prop::Symbol; kwargs...)

Sort the particles in-place according to the specified property.

All property arrays of the particles are rearranged in-place according to the property being sorted.
The keyword arguments are passed on to [`Base.sortperm`](https://docs.julialang.org/en/v1/base/sort/#Base.sortperm).

The sorting algorithm for unitful properties may also be `RadixSort` from
[SortingAlgorithms.jl](https://github.com/JuliaCollections/SortingAlgorithms.jl).
"""
function Base.sort!(p::AbstractParticles, prop::Symbol; kwargs...)
    vals = p[prop]
    vals isa Number && return p
    ind = usortperm(vals; kwargs...)
    return applyind!(p, ind)
end

"""
    Base.sort!(pc::AbstractParticleCollection, prop::Symbol; kwargs...)

Sort the particles in the collection in-place by calling `Base.sort!` on each of the [`Particles`](@ref) objects.
"""
function Base.sort!(pc::AbstractParticleCollection, prop::Symbol; kwargs...)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        sort!(pc[ptype], prop; kwargs...)
    end

    return pc
end

"""
    Base.sort(p::AbstractParticles, prop::Symbol; affect=keys(p), kwargs...)

Create new particles with the particles sorted according to the specified property.

All property arrays of the particles are rearranged according to the property being sorted.
If the keyword argument `affect` is a non-empty vector of `Symbol`s, only those properties are rearranged and added
to the newly created particles object.
The other keyword arguments are passed on to [`Base.sortperm`](https://docs.julialang.org/en/v1/base/sort/#Base.sortperm).

The sorting algorithm for unitful properties may also be `RadixSort` from
[SortingAlgorithms.jl](https://github.com/JuliaCollections/SortingAlgorithms.jl).
"""
function Base.sort(p::AbstractParticles, prop::Symbol; affect=keys(p), kwargs...)
    vals = p[prop]
    vals isa Number && return applyind(p, (:); affect)
    ind = usortperm(vals; kwargs...)
    return applyind(p, ind; affect)
end

"""
    Base.sort(pc::AbstractParticleCollection, prop::Symbol; [affect,] kwargs...)

Create new particle collection with the particles in the collection sorted by calling `Base.sort` on each of the particles objects.

If the keyword argument `affect` is a non-empty vector of `Symbol`s, only the affected properties are kept for the particles.
The specified affected properties do not have to be available for all particles.
For collections of [`Particles`](@ref), `affect` can alternatively be a vector of tuples or pairs in the following form:
`[(:dm, [:id, :pos, :mass]), (:gas, [:id, :pos, :mass, :temp])]`,
`[:dm => [:id, :pos, :mass], :gas => [:id, :pos, :mass, :temp]]`.
"""
function Base.sort(pc::AbstractParticleCollection, prop::Symbol; affect=nothing, kwargs...)
    # preallocate particle collection (to avoid threading issues)
    pcnew = _preallocate_particle_collection(pc, affect)

    if affect isa AbstractVector{<:Union{Tuple,Pair}}
        Threads.@threads for (ptype, props) in affect |> collect # collect used for compatibility with threads
            pcnew[ptype] = sort(pc[ptype], prop; affect=props, kwargs...)
        end
    else
        Threads.@threads for ptype in keys(pc) |> collect
            if isnothing(affect)
                pcnew[ptype] = sort(pc[ptype], prop; kwargs...)
            else
                pcnew[ptype] = sort(pc[ptype], prop; affect, kwargs...)
            end
        end
    end

    return pcnew
end


"""
    Base.filter!(f, p::AbstractParticles)
    Base.filter!(f, pc::AbstractParticleCollection)

Filter the particles or collection in-place by masks returned by the function `f`.

The function `f` takes a particles object as argument and has to either return a `BitArray` mask with a length of the
number of particles or an array of indices.
"""
function Base.filter!(f, p::AbstractParticles)
    ind = f(p)
    return applyind!(p, ind)
end

function Base.filter!(f, pc::AbstractParticleCollection)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        filter!(f, pc[ptype])
    end

    return pc
end

"""
    Base.filter(f, p::AbstractParticles; affect=keys(p))
    Base.filter(f, pc::AbstractParticleCollection[; affect])

Create new particles or collection with the particles filtered by a mask returned by the function `f`.

The function `f` takes the particles as argument and has to either return a `BitArray` mask with a length of the
number of particles or an array of indices.
If the keyword argument `affect` is a non-empty vector of `Symbol`s, only those properties are filtered and added
to the newly created particles object.

For collections of [`Particles`](@ref), `affect` can alternatively be a vector of tuples or pairs in the following form:
`[(:dm, [:id, :pos, :mass]), (:gas, [:id, :pos, :mass, :temp])]`,
`[:dm => [:id, :pos, :mass], :gas => [:id, :pos, :mass, :temp]]`.
"""
function Base.filter(f, p::AbstractParticles; affect=keys(p))
    ind = f(p)
    return applyind(p, ind; affect)
end

function Base.filter(f, pc::AbstractParticleCollection; affect=nothing)
    # preallocate particle collection (to avoid threading issues)
    pcnew = _preallocate_particle_collection(pc, affect)

    if affect isa AbstractVector{<:Union{Tuple,Pair}}
        Threads.@threads for (ptype, props) in affect |> collect # collect used for compatibility with threads
            pcnew[ptype] = filter(f, pc[ptype]; affect=props)
        end
    else
        Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
            if isnothing(affect)
                pcnew[ptype] = filter(f, pc[ptype])
            else
                pcnew[ptype] = filter(f, pc[ptype]; affect)
            end
        end
    end

    return pcnew
end


"""
    Base.filter!(p::AbstractParticles; ids, idprop=:id)
    Base.filter!(pc::AbstractParticleCollection; ids, idprop=:id)

Filter the particles or collection in-place by keeping only the particles with the given IDs.
"""
function Base.filter!(p::AbstractParticles; ids, idprop=:id)
    ind = findall_in(p[idprop], ids)
    return applyind!(p, ind)
end

function Base.filter!(pc::AbstractParticleCollection; ids, idprop=:id)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        filter!(pc[ptype]; ids, idprop)
    end

    return pc
end

"""
    Base.filter(p::AbstractParticles; ids, affect=keys(p), idprop=:id)
    Base.filter(pc::AbstractParticleCollection; ids, idprop=:id[, affect])

Create new particles or collection with them filtered by keeping only the particles with the given IDs.

If the keyword argument `affect` is a non-empty vector of `Symbol`s, only those properties are filtered and added
to the newly created particles object.

For collections of [`Particles`](@ref), `affect` can alternatively be a vector of tuples or pairs in the following form:
`[(:dm, [:id, :pos, :mass]), (:gas, [:id, :pos, :mass, :temp])]`,
`[:dm => [:id, :pos, :mass], :gas => [:id, :pos, :mass, :temp]]`.
"""
function Base.filter(p::AbstractParticles; ids, idprop=:id, affect=keys(p))
    ind = findall_in(p[idprop], ids)
    return applyind(p, ind; affect)
end

function Base.filter(pc::AbstractParticleCollection; ids, idprop=:id, affect=nothing)
    # preallocate particle collection (to avoid threading issues)
    pcnew = _preallocate_particle_collection(pc, affect)

    if affect isa AbstractVector{<:Union{Tuple,Pair}}
        Threads.@threads for (ptype, props) in affect |> collect # collect used for compatibility with threads
            pcnew[ptype] = filter(pc[ptype]; ids, idprop, affect=props)
        end
    else
        Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
            if isnothing(affect)
                pcnew[ptype] = filter(pc[ptype]; ids, idprop)
            else
                pcnew[ptype] = filter(pc[ptype]; ids, idprop, affect)
            end
        end
    end

    return pcnew
end


"""
    Base.filter(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos; affect=keys(p))
    Base.filter(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos[; affect])

Create new particles or collection with them filtered by keeping only the particles inside the given [`AbstractCosmoGeometry`](@ref).

The geometry filter is applied to the property specified.
If the keyword argument `affect` is a non-empty vector of `Symbol`s, only those properties are filtered and added
to the newly created particles object.

For collections of [`Particles`](@ref), `affect` can alternatively be a vector of tuples or pairs in the following form:
`[(:dm, [:id, :pos, :mass]), (:gas, [:id, :pos, :mass, :temp])]`,
`[:dm => [:id, :pos, :mass], :gas => [:id, :pos, :mass, :temp]]`.
"""
function Base.filter(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos; affect=keys(p))
    ind = mask_in(p[prop], geo)
    return applyind(p, ind; affect)
end

function Base.filter(
    pc::AbstractParticleCollection,
    geo::AbstractCosmoGeometry,
    prop::Symbol=:pos;
    affect=nothing,
)
    # preallocate particle collection (to avoid threading issues)
    pcnew = _preallocate_particle_collection(pc, affect)

    if affect isa AbstractVector{<:Union{Tuple,Pair}}
        Threads.@threads for (ptype, props) in affect |> collect # collect used for compatibility with threads
            pcnew[ptype] = filter(pc[ptype], geo, prop; affect=props)
        end
    else
        Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
            if isnothing(affect)
                pcnew[ptype] = filter(pc[ptype], geo, prop)
            else
                pcnew[ptype] = filter(pc[ptype], geo, prop; affect)
            end
        end
    end

    return pcnew
end

"""
    Base.filter!(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    Base.filter!(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos)

Filter the particles or collection in-place by keeping only the particles inside the given [`AbstractCosmoGeometry`](@ref).

The geometry filter is applied to the property specified.
"""
function Base.filter!(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    ind = mask_in(p[prop], geo)
    return applyind!(p, ind)
end

function Base.filter!(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        filter!(pc[ptype], geo, prop)
    end

    return pc
end


"""
    Base.delete!(p::AbstractParticles; ids)
    Base.delete!(pc::AbstractParticleCollection; ids)

Filter the particles or collection in-place by removing all particles with the given IDs.
"""
function Base.delete!(p::AbstractParticles; ids)
    ind = findall_in(p.id, ids)
    return removeind!(p, ind)
end

function Base.delete!(pc::AbstractParticleCollection; ids)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        delete!(pc[ptype]; ids)
    end

    return pc
end

"""
    delete(p::AbstractParticles, ids; affect=keys(p))
    delete(pc::AbstractParticleCollection; ids[, affect])

Create new particles or collection with them filtered by keeping only the particles with the given IDs.

If the keyword argument `affect` is a non-empty vector of `Symbol`s, only those properties are filtered and added
to the newly created particles object.

For collections of [`Particles`](@ref), `affect` can alternatively be a vector of tuples or pairs in the following form:
`[(:dm, [:id, :pos, :mass]), (:gas, [:id, :pos, :mass, :temp])]`,
`[:dm => [:id, :pos, :mass], :gas => [:id, :pos, :mass, :temp]]`.
"""
function delete(p::AbstractParticles; ids, affect=keys(p))
    ind = findall_in(p.id, ids)
    return removeind(p, ind; affect)
end

function delete(pc::AbstractParticleCollection; ids, affect=nothing)
    # preallocate particle collection (to avoid threading issues)
    pcnew = _preallocate_particle_collection(pc, affect)

    if affect isa AbstractVector{<:Union{Tuple,Pair}}
        Threads.@threads for (ptype, props) in affect |> collect # collect used for compatibility with threads
            pcnew[ptype] = delete(pc[ptype]; ids, affect=props)
        end
    else
        Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
            if isnothing(affect)
                pcnew[ptype] = delete(pc[ptype]; ids)
            else
                pcnew[ptype] = delete(pc[ptype]; ids, affect)
            end
        end
    end

    return pcnew
end


"""
    Base.delete(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos; affect=keys(p))
    Base.delete(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos[; affect])

Create new particles or collection with them filtered by keeping only the particles outside the given [`AbstractCosmoGeometry`](@ref).

The geometry filter is applied to the property specified.
If the keyword argument `affect` is a non-empty vector of `Symbol`s, only those properties are filtered and added
to the newly created particles object.

For collections of [`Particles`](@ref), `affect` can alternatively be a vector of tuples or pairs in the following form:
`[(:dm, [:id, :pos, :mass]), (:gas, [:id, :pos, :mass, :temp])]`,
`[:dm => [:id, :pos, :mass], :gas => [:id, :pos, :mass, :temp]]`.
"""
function delete(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos; affect=keys(p))
    ind = .~mask_in(p[prop], geo)
    return applyind(p, ind; affect)
end

function delete(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos; affect=nothing)
    # preallocate particle collection (to avoid threading issues)
    pcnew = _preallocate_particle_collection(pc, affect)

    if affect isa AbstractVector{<:Union{Tuple,Pair}}
        Threads.@threads for (ptype, props) in affect |> collect # collect used for compatibility with threads
            pcnew[ptype] = delete(pc[ptype], geo, prop; affect=props)
        end
    else
        Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
            if isnothing(affect)
                pcnew[ptype] = delete(pc[ptype], geo, prop)
            else
                pcnew[ptype] = delete(pc[ptype], geo, prop; affect)
            end
        end
    end

    return pcnew
end

"""
    Base.delete!(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    Base.delete!(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos)

Filter the particles or collection in-place by keeping only the particles outside the given [`AbstractCosmoGeometry`](@ref).

The geometry filter is applied to the property specified.
"""
function Base.delete!(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    ind = .~mask_in(p[prop], geo)
    return applyind!(p, ind)
end

function Base.delete!(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    Threads.@threads for ptype in keys(pc) |> collect # collect used for compatibility with threads
        delete!(pc[ptype], geo, prop)
    end

    return pc
end



function _preallocate_particle_collection(pc, affect)
    pcnew = empty(pc)
    if affect isa AbstractVector{<:Union{Tuple,Pair}}
        for (ptype, _) in affect
            pcnew[ptype] = pc[ptype]
        end
    else
        for ptype in keys(pc)
            pcnew[ptype] = pc[ptype]
        end
    end
    return pcnew
end
