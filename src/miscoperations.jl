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
    Base.sort(p::AbstractParticles, prop::Symbol; affect=(), kwargs...)

Create new particles with the particles sorted according to the specified property.

All property arrays of the particles are rearranged according to the property being sorted.
If the keyword argument `affect` is a non-empty tuple of `Symbol`s, only those properties are rearranged and added
to the newly created particles object (besides the originally sorted property).
The other keyword arguments are passed on to [`Base.sortperm`](https://docs.julialang.org/en/v1/base/sort/#Base.sortperm).

The sorting algorithm for unitful properties may also be `RadixSort` from
[SortingAlgorithms.jl](https://github.com/JuliaCollections/SortingAlgorithms.jl).
"""
function Base.sort(p::AbstractParticles, prop::Symbol; affect=(), kwargs...)
    ind = usortperm(p[prop]; kwargs...)
    return applyind(p, ind; affect=ifelse(isempty(affect), (), union(affect, (prop,))))
end


"""
    Base.filter!(f, p::AbstractParticles)

Filter the particles in-place by a mask returned by the function `f`.

The function `f` takes the particles as argument and has to either return a `BitArray` mask with a length of the
number of particles or an array of indices.
"""
function Base.filter!(f, p::AbstractParticles)
    ind = f(p)
    return applyind!(p, ind)
end

"""
    Base.filter(f, p::AbstractParticles; affect=())

Create new particles with them filtered by a mask returned by the function `f`.

The function `f` takes the particles as argument and has to either return a `BitArray` mask with a length of the
number of particles or an array of indices.
If the keyword argument `affect` is a non-empty tuple of `Symbol`s, only those properties are filtered and added
to the newly created particles object.
"""
function Base.filter(f, p::AbstractParticles; affect=())
    ind = f(p)
    return applyind(p, ind; affect)
end


"""
    Base.filter!(p::AbstractParticles, ids)

Filter the particles in-place by keeping only those with the given IDs.
"""
function Base.filter!(p::AbstractParticles; ids)
    ind = findall_in(p.id, ids)
    return applyind!(p, ind)
end

"""
    Base.filter(p::AbstractParticles, ids; affect=())

Create new particles with them filtered by keeping only those with the given IDs.

If the keyword argument `affect` is a non-empty tuple of `Symbol`s, only those properties are filtered and added
to the newly created particles object.
"""
function Base.filter(p::AbstractParticles; ids, affect=())
    ind = findall_in(p.id, ids)
    return applyind(p, ind; affect)
end
