"""
    CosmoParticles._applyind(a, ind::AbstractVector)

Apply indices or mask to a `Number`, `Vector`, or `Matrix`.

The following indexing is applied:
- `a::Number`: `a` is returned directly.
- `a::Vector`: `a[ind]` is returned.
- `a::Matrix`: `a[:, ind]` is returned.

This is not exported.
"""
_applyind(a::Number, _::AbstractVector) = a
_applyind(a::AbstractVector, ind::AbstractVector) = a[ind]
_applyind(a::AbstractMatrix, ind::AbstractVector) = a[:, ind]

"""
    CosmoParticles.applyind!(p::AbstractParticles, ind::AbstractVector)

In-place application of indices or mask to all particle properties.

This is not exported.
"""
function applyind!(p::AbstractParticles, ind::AbstractVector)
    for key in keys(p)
        p[key] = _applyind(p[key], ind)
    end

    return p
end

"""
    CosmoParticles.applyind(p::AbstractParticles, ind::AbstractVector; affect=keys(p))

Create new particles with the given indices or mask applied to all particle properties.

This can also be called by the simple syntax `p[ind]`.
If the keyword argument `affect` is a non-empty tuple of `Symbol`s, only those properties are indexed into and
added to the newly created particles object.

This is not exported.
"""
function applyind(p::AbstractParticles, ind::AbstractVector; affect=keys(p))
    pnew = empty(p)

    for key in intersect(affect, keys(p))
        pnew[key] = _applyind(p[key], ind)
    end

    return pnew
end

"""
    CosmoParticles.findall_in(a::AbstractVector, set)

Return all indices of `a` that are in `set`.

If both `a` and `set` are sorted `AbstractVector`s, then the optimized [`findall_in_sorted`](@ref) is called.
Otherwise, a `Set` is constructed from the `Vector` to perform the checks with `in`.

This is not exported.
"""
function findall_in(a::AbstractVector, set::AbstractVector)
    if issorted(a) && issorted(set)
        return findall_in_sorted(a, set)
    else
        return findall_in(a, Set(set))
    end
end

findall_in(a::AbstractVector, set::AbstractSet) = findall(in.(a, (set,)))

"""
    CosmoParticles.findall_in_sorted(a::AbstractVector, set::AbstractVector)

Return all indices of `a` that are in `set`, where both `a` and `set` are assumed to be sorted.

This uses an optimized algorithm that is faster than creating a `Set` from `set` and performing checks with `in`.

This is not exported.
"""
function findall_in_sorted(a::AbstractVector, set::AbstractVector)
    if isempty(a) || isempty(set)
        return Int64[]
    end

    na = length(a)
    nset = length(set)

    ind_all = Vector{Int64}(undef, min(na, nset))

    # NOTE: this is an alternative algorithm, which performs slightly worse than the below
    # iind = 0
    # iset = 1
    # for ia in 1:na
    # while a[ia] > set[iset] && iset < nset
    # iset += 1
    # end
    # if a[ia] == set[iset]
    # iind += 1
    # iset += 1
    # ind_all[iind] = ia
    # end
    # end

    iind = 0
    ia = 1
    iset = 1

    while true
        if a[ia] == set[iset]
            iind += 1
            ind_all[iind] = ia
            iset += 1
            ia += 1
            ((ia > na) || (iset > nset)) && break
        else
            if a[ia] < set[iset]
                ia += 1
                ia > na && break
            else
                iset += 1
                iset > nset && break
            end
        end
    end

    return resize!(ind_all, iind)
end

"""
    product_preserve_type(arr::AbstractArray{T}, b::Real) where {T}

Multiply an array with a scalar while preserving the element type of the array.

For unitful arrays, the scalar factor is converted to the number type of the quantity before multiplying.
"""
product_preserve_type(arr::AbstractArray{T}, b::Real) where {T<:Real} = arr .* convert(T, b)
product_preserve_type(arr::AbstractArray{Quantity{T}}, b::Real) where {T<:Real} = arr .* convert(T, b)

"""
    product_preserve_type!(arr::AbstractArray{T}, b::Real) where {T}

Multiply an array with a scalar in-place while preserving the element type of the array.

By converting the scalar factor to the array element type, this can be more performant than a normal
broadcasted product.
For unitful arrays, the scalar factor is converted to the number type of the quantity before multiplying.
"""
function product_preserve_type!(arr::AbstractArray{T}, b::Real) where {T<:Real}
    arr .*= convert(T, b)
end
function product_preserve_type!(arr::AbstractArray{Quantity{T}}, b::Real) where {T<:Real}
    arr .*= convert(T, b)
end
