"""
    Base.sort!(p::AbstractParticles, prop::Symbol; kwargs...)

Sort the particles in-place according to the specified property.

All property arrays of the particles are rearranged in-place according to the property being sorted.
The keyword arguments are passed on to [`Base.sortperm`](https://docs.julialang.org/en/v1/base/sort/#Base.sortperm).

The sorting algorithm for unitful properties may also be `RadixSort` from
[SortingAlgorithms.jl](https://github.com/JuliaCollections/SortingAlgorithms.jl).
"""
function Base.sort!(p::AbstractParticles, prop::Symbol; kwargs...)
    ind = usortperm(p[prop]; kwargs...)
    return applyind!(p, ind)
end

"""
    Base.sort!(pc::AbstractParticleCollection, prop::Symbol; kwargs...)

Sort the particles in the collection in-place by calling `Base.sort!` on each of the [`Particles`](@ref) objects.
"""
function Base.sort!(pc::AbstractParticleCollection, prop::Symbol; kwargs...)
    for ptype in keys(pc)
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
    ind = usortperm(p[prop]; kwargs...)
    return applyind(p, ind; affect)
end

"""
    Base.sort(pc::AbstractParticleCollection, prop::Symbol; [affect,] kwargs...)

Create new particle collection with the particles in the collection sorted by calling `Base.sort` on each of the [`Particles`](@ref) objects.

If the keyword argument `affect` is a non-empty vector of `Symbol`s, only the affected properties are kept for the particles.
The specified affected properties do not have to be available for all particles.
"""
function Base.sort(pc::AbstractParticleCollection, prop::Symbol; affect=nothing, kwargs...)
    pc = copy(pc)

    for ptype in keys(pc)
        if isnothing(affect)
            pc[ptype] = sort(pc[ptype], prop; kwargs...)
        else
            pc[ptype] = sort(pc[ptype], prop; affect, kwargs...)
        end
    end

    return pc
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
    for ptype in keys(pc)
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
"""
function Base.filter(f, p::AbstractParticles; affect=keys(p))
    ind = f(p)
    return applyind(p, ind; affect)
end

function Base.filter(f, pc::AbstractParticleCollection; affect=nothing)
    pc = copy(pc)

    for ptype in keys(pc)
        if isnothing(affect)
            pc[ptype] = filter(f, pc[ptype])
        else
            pc[ptype] = filter(f, pc[ptype]; affect)
        end
    end

    return pc
end


"""
    Base.filter!(p::AbstractParticles; ids)
    Base.filter!(pc::AbstractParticleCollection; ids)

Filter the particles or collection in-place by keeping only the particles with the given IDs.
"""
function Base.filter!(p::AbstractParticles; ids)
    ind = findall_in(p.id, ids)
    return applyind!(p, ind)
end

function Base.filter!(pc::AbstractParticleCollection; ids)
    for ptype in keys(pc)
        filter!(pc[ptype]; ids)
    end

    return pc
end

"""
    Base.filter(p::AbstractParticles, ids; affect=keys(p))
    Base.filter(pc::AbstractParticleCollection; ids[, affect])

Create new particles or collection with them filtered by keeping only the particles with the given IDs.

If the keyword argument `affect` is a non-empty vector of `Symbol`s, only those properties are filtered and added
to the newly created particles object.
"""
function Base.filter(p::AbstractParticles; ids, affect=keys(p))
    ind = findall_in(p.id, ids)
    return applyind(p, ind; affect)
end

function Base.filter(pc::AbstractParticleCollection; ids, affect=nothing)
    pc = copy(pc)

    for ptype in keys(pc)
        if isnothing(affect)
            pc[ptype] = filter(pc[ptype]; ids)
        else
            pc[ptype] = filter(pc[ptype]; ids, affect)
        end
    end

    return pc
end


"""
    Base.filter(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos; affect=keys(p))
    Base.filter(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos[; affect])

Create new particles or collection with them filtered by keeping only the particles inside the given [`AbstractCosmoGeometry`](@ref).

The filter is applied to the property specified.
If the keyword argument `affect` is a non-empty vector of `Symbol`s, only those properties are filtered and added
to the newly created particles object.
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
    pc = copy(pc)

    for ptype in keys(pc)
        if isnothing(affect)
            pc[ptype] = filter(pc[ptype], geo, prop)
        else
            pc[ptype] = filter(pc[ptype], geo, prop; affect)
        end
    end

    return pc
end

"""
    Base.filter!(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    Base.filter!(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos)

Filter the particles or collection in-place by keeping only the particles inside the given [`AbstractCosmoGeometry`](@ref).

The filter is applied to the property specified.
"""
function Base.filter!(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    ind = mask_in(p[prop], geo)
    return applyind!(p, ind)
end

function Base.filter!(pc::AbstractParticleCollection, geo::AbstractCosmoGeometry, prop::Symbol=:pos)
    for ptype in keys(pc)
        filter!(pc[ptype], geo, prop)
    end

    return pc
end
