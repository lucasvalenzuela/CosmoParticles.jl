"""
    get_dims(projection::Symbol)
    get_dims(dims::AbstractArray)

Returns the 2D dimensions of the passed projection, or the dimensions themselves.

The following dimensions are returned:
- `[1, 3]`: `:edgeon` or `:xz`
- `[2, 3]`: `:sideon` or `:yz`
- `[1, 2]`: `:faceon` or `:xy`
"""
function get_dims(projection::Symbol)
    if projection === :edgeon || projection === :xz
        return [1, 3]
    elseif projection === :sideon || projection === :yz
        return [2, 3]
    elseif projection === :faceon || projection === :xy
        return [1, 2]
    else
        @error "The projection `:$projection` is not valid. Please use `:edgeon`, `:sideon`, `:faceon`, `:xy`, `:yz`, or `:xz`."
    end
end
get_dims(dims) = dims

"""
    get_losdim(projection::Symbol)
    get_losdim(dims::AbstractArray)

Returns the 1D line-of-sight dimension of the passed projection, or of the passed 2D dimensions.

The following dimensions are returned:
- `2`: `:edgeon` or `:xz` or `[1, 3]`
- `1`: `:sideon` or `:yz` or `[2, 3]`
- `3`: `:faceon` or `:xy` or `[1, 2]`
"""
function get_losdim(projection::Symbol)
    if projection === :edgeon || projection === :xz
        return 2
    elseif projection === :sideon || projection === :yz
        return 1
    elseif projection === :faceon || projection === :xy
        return 3
    else
        @error "The projection `:$projection` is not valid. Please use `:edgeon`, `:sideon`, `:faceon`, `:xy`, `:yz`, or `:xz`."
    end
end
get_losdim(dims) = setdiff(1:3, dims)[1]



"""
    CosmoParticles._applyind(a, ind::Union{AbstractVector,Colon})

Apply indices or mask to a `Number`, `Vector`, or `Matrix`.

Other possible types: `Nothing`

The following indexing is applied:
- `a::Number`: `a` is returned directly.
- `a::Vector`: `a[ind]` is returned.
- `a::Matrix`: `a[:, ind]` is returned.

This is not exported.
"""
_applyind(a::Union{Number,Nothing}, _::Union{AbstractVector,Colon}) = a
_applyind(a::AbstractVector, ind::Union{AbstractVector,Colon}) = a[ind]
_applyind(a::AbstractMatrix, ind::Union{AbstractVector,Colon}) = a[:, ind]
_applyind(a::ApplyVector{<:Any,F}, mask::AbstractVector{Bool}) where {F<:Union{typeof(hcat),typeof(vcat)}} =
    a[mask]
_applyind(a::ApplyMatrix{<:Any,F}, mask::AbstractVector{Bool}) where {F<:Union{typeof(hcat),typeof(vcat)}} =
    a[:, mask]

function _applyind(a::ApplyVector{<:Any,F}, ind::AbstractVector) where {F<:Union{typeof(hcat),typeof(vcat)}}
    anew = similar(a, length(ind))
    @inbounds for (i, indi) in enumerate(ind)
        anew[i] = a[indi]
    end
    return anew
end

function _applyind(a::ApplyMatrix{<:Any,F}, ind::AbstractVector) where {F<:Union{typeof(hcat),typeof(vcat)}}
    n = size(a, 1)
    anew = similar(a, n, length(ind))
    @inbounds for (i, indi) in enumerate(ind)
        for j in 1:n
            anew[j, i] = a[j, indi]
        end
    end
    return anew
end

"""
    CosmoParticles.applyind!(p::AbstractParticles, ind::Union{AbstractVector,Colon})

In-place application of indices or mask to all particle properties.

This is not exported.
"""
function applyind!(p::AbstractParticles, ind::Union{AbstractVector,Colon})
    Threads.@threads for key in keys(p) |> collect # collect used for compatibility with threads
        p[key] = _applyind(p[key], ind)
    end

    return p
end

"""
    CosmoParticles.applyind(p::AbstractParticles, ind::Union{AbstractVector,Colon}; affect=keys(p))

Create new particles with the given indices or mask applied to all particle properties.

This can also be called by the simple syntax `p[ind]`.
If the keyword argument `affect` is a vector of `Symbol`s, only those properties are indexed into and
added to the newly created particles object.

This is not exported.
"""
function applyind(p::AbstractParticles, ind::Union{AbstractVector,Colon}; affect=keys(p))
    affected_keys = intersect(affect, keys(p)) |> collect # collect used for compatibility with threads

    # preallocate particles (to avoid threading issues)
    pnew = empty(p)
    for key in affected_keys
        pnew[key] = p[key]
    end

    Threads.@threads for key in affected_keys # collect used for compatibility with threads
        pnew[key] = _applyind(p[key], ind)
    end

    return pnew
end




"""
    CosmoParticles._removeind(a, ind::AbstractVector)

Remove indices of a `Number`, `Vector`, or `Matrix`.

The vector `ind` has to contain sorted and unique indices: `sort!(unique!(ind))`.
Additionally, all elements in `ind` have to be valid indices in `a`.

Other possible types: `Nothing`

The following indexing is applied:
- `a::Number`: `a` is returned directly.
- `a::Vector`: `a[Not(ind)]` is returned.
- `a::Matrix`: `a[:, Not(ind)]` is returned.

This is not exported.
"""
_removeind(a::Union{Number,Nothing}, _::AbstractVector) = a

function _removeind(a::AbstractVector, ind::AbstractVector)
    na = length(a)
    nind = length(ind)
    nanew = na - nind
    anew = similar(a, nanew)

    # catch empty array (otherwise throws BoundsError in loop)
    na == nind && return anew

    ia = 1
    iind = 1
    ianew = 1
    @inbounds while true
        while iind > nind || ia < ind[iind]
            anew[ianew] = a[ia]
            ia += 1
            ianew += 1
            ianew > nanew && return anew
        end
        ia += 1
        iind += 1
    end
end

function _removeind(a::AbstractMatrix, ind::AbstractVector)
    n = size(a, 1)
    na = size(a, 2)
    nind = length(ind)
    nanew = na - nind
    anew = similar(a, n, nanew)

    # catch empty array (otherwise throws BoundsError in loop)
    na == nind && return anew

    ia = 1
    iind = 1
    ianew = 1
    @inbounds while true
        while iind > nind || ia < ind[iind]
            for j in 1:n
                anew[j, ianew] = a[j, ia]
            end
            ia += 1
            ianew += 1
            ianew > nanew && return anew
        end
        ia += 1
        iind += 1
    end
end

"""
    CosmoParticles.removeind!(p::AbstractParticles, ind::AbstractVector)

In-place removal of indices to all particle properties.

This is not exported.
"""
function removeind!(p::AbstractParticles, ind::AbstractVector)
    Threads.@threads for key in keys(p) |> collect # collect used for compatibility with threads
        p[key] = _removeind(p[key], ind)
    end

    return p
end

"""
    CosmoParticles.removeind(p::AbstractParticles, ind::AbstractVector; affect=keys(p))

Create new particles without the given indices or mask applied to all particle properties.

This can also be called by the simple syntax `p[ind]`.
If the keyword argument `affect` is a vector of `Symbol`s, only those properties are indexed into and
added to the newly created particles object.

This is not exported.
"""
function removeind(p::AbstractParticles, ind::AbstractVector; affect=keys(p))
    affected_keys = intersect(affect, keys(p)) |> collect # collect used for compatibility with threads

    # preallocate particles (to avoid threading issues)
    pnew = empty(p)
    for key in affected_keys
        pnew[key] = p[key]
    end

    Threads.@threads for key in affected_keys
        pnew[key] = _removeind(p[key], ind)
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
    iset = 1
    ia = 1

    # for a large number of values, search for the first indices in a faster way
    if na > 1_000_000
        if first(a) < first(set)
            ia = searchsortedfirst(a, first(set))
        else
            iset = searchsortedfirst(set, first(a))
            ia = searchsortedfirst(a, set[iset])
        end
    end

    if ia > na || iset > nset
        return Int64[]
    end

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
    CosmoParticles.product_preserve_type(arr::AbstractArray{T}, b::Real) where {T}

Multiply an array with a scalar while preserving the element type of the array.

For unitful arrays, the scalar factor is converted to the number type of the quantity before multiplying.

This is not exported.
"""
product_preserve_type(arr::AbstractArray{T}, b::Real) where {T<:Real} = arr .* convert(T, b)
product_preserve_type(arr::AbstractArray{<:Quantity{T}}, b::Real) where {T<:Real} = arr .* convert(T, b)

"""
    CosmoParticles product_preserve_type!(arr::AbstractArray{T}, b::Real) where {T}

Multiply an array with a scalar in-place while preserving the element type of the array.

By converting the scalar factor to the array element type, this can be more performant than a normal
broadcasted product.
For unitful arrays, the scalar factor is converted to the number type of the quantity before multiplying.

This is not exported.
"""
function product_preserve_type!(arr::AbstractArray{T}, b::Real) where {T<:Real}
    arr .*= convert(T, b)
end
function product_preserve_type!(arr::AbstractArray{<:Quantity{T}}, b::Real) where {T<:Real}
    arr .*= convert(T, b)
end
function product_preserve_type!(
    arr::ApplyArray{<:Any,N,F},
    b::Real,
) where {N,F<:Union{typeof(hcat),typeof(vcat)}}
    product_preserve_type!.(arr.args, b)
    return arr
end


"""
    CosmoParticles.ustrip_lazy([unit,] a::Number)
    CosmoParticles.ustrip_lazy([unit,] a::AbstractArray)
    CosmoParticles.ustrip_lazy!(unit, a::AbstractArray)

Strips off the units of unitful arrays in a performant way while keeping the number type conserved (not always the case for `ustrip`).

Like [`ustrip`](https://painterqubits.github.io/Unitful.jl/stable/manipulations/#Unitful.ustrip) for normal
unitful arrays, but without reallocating non-unitful or lazy arrays. Differently to `ustrip`, this also works
for `Number`s and arrays thereof, simply returning the original array.

Note that after calling this method with a different conversion unit, the original array still has the old units,
but the numerical values correspond to the new units. This method should always be called as
`a = uconvert_lazy!(u, a)`.

This is not exported.
"""
ustrip_lazy(a::Unitful.AbstractQuantity) = ustrip(a)
ustrip_lazy(a::Number) = a
ustrip_lazy(a::AbstractArray) = a
ustrip_lazy(a::Array{<:Quantity}) = ustrip(a)
ustrip_lazy(a::Diagonal{<:Quantity}) = ustrip(a)
ustrip_lazy(a::Bidiagonal{<:Quantity}) = ustrip(a)
ustrip_lazy(a::Tridiagonal{<:Quantity}) = ustrip(a)
ustrip_lazy(a::SymTridiagonal{<:Quantity}) = ustrip(a)

ustrip_lazy(a::AbstractArray{<:Quantity{T}}) where {T} = reinterpret(T, a) #ustrip.(a) # alternative implementation
ustrip_lazy(a::ApplyArray{<:Quantity{T}}) where {T} = reinterpret(T, a)
ustrip_lazy(a::Fill{<:Quantity{T}}) where {T} = reinterpret(T, a)

ustrip_lazy!(u::Unitful.Units, a) = ustrip_lazy(uconvert_lazy!(u, a))
ustrip_lazy(_::Unitful.Units, a::Number) = a
ustrip_lazy(_::Unitful.Units, a::AbstractArray{<:Number}) = a
ustrip_lazy(u::Unitful.Units, a::AbstractArray{<:Quantity}) = ustrip.(u, a)
# ustrip_lazy(u::Unitful.Units, a::AbstractVector{<:Quantity{T}}) where {T<:AbstractFloat} = ustrip_lazy.(u, a)
ustrip_lazy(u::Unitful.Units, a::Quantity{T}) where {T<:Real} = ustrip_lazy(uconvert_lazy(u, a))
function ustrip_lazy(u::Unitful.Units, a::AbstractArray{<:Quantity{T}}) where {T<:AbstractFloat} 
    unit(eltype(a)) == u ? ustrip_lazy(a) : ustrip_lazy.(u, a)
end

uconvert_lazy(u::Unitful.Units, a::AbstractArray{<:Quantity}) = uconvert_lazy!(u, copy(a))

"""
    CosmoParticles.uconvert_lazy(u::Unitful.Units, a::AbstractArray{<:Quantity})

Converts the unitful array to the given units while preserving the number type passed (not always the case for `uconvert`).

This is not exported.
"""
function uconvert_lazy(u::Unitful.Units, a::Quantity{T,D,U}) where {T,D,U}
    return convert(T, ustrip(u, a)) * u
end

"""
    CosmoParticles.uconvert_lazy!(u::Unitful.Units, a::AbstractArray{<:Quantity{T,D,U}}) where {T,D,U}

Converts the unitful array to the specified units in-place.

Note that after calling this method the original array still has the old units, but the numerical values
correspond to the new units. This method should always be called as `a = uconvert_lazy!(u, a)`.

This is not exported.
"""
function uconvert_lazy!(u::Unitful.Units, a::AbstractArray{<:Quantity{T,D,U}}) where {T,D,U}
    Q = Quantity{T,dimension(u),typeof(u)}
    if typeof(u) == U
        return reinterpret(Q, a)
    elseif (u isa Unitful.AffineUnits) || (U <: Unitful.AffineUnits)
        a_unitless = ustrip_lazy(a)
        a_unitless .= ustrip_lazy(Unitful.uconvert_affine.(u, a))
        return reinterpret(Q, a_unitless)
    else
        product_preserve_type!(a, Unitful.convfact(u, U()))
        return reinterpret(Q, a)
    end
end
