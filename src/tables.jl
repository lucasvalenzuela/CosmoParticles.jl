Tables.istable(::Type{<:AbstractParticles}) = true

function Tables.schema(p::AbstractParticles) 
    t = _get_eltype.(values(p))
    return Tables.Schema(Tables.columnnames(p), t)
end

_get_eltype(::T) where {T} = T
_get_eltype(::AbstractVector{T}) where {T} = T
_get_eltype(::AbstractMatrix{T}) where {T} = Vector{T}

Tables.columnaccess(::Type{<:AbstractParticles}) = true
Tables.columns(p::AbstractParticles) = p

function Tables.getcolumn(p::AbstractParticles, ::Type{T}, col::Int, nm::Symbol) where {T}
    val = p[nm]
    return val isa AbstractVector ? val : fill(val, nrow(p))
end

function Tables.getcolumn(p::AbstractParticles, ::Type{<:Vector}, col::Int, nm::Symbol)
    eachcol(p[nm])
end

function Tables.getcolumn(p::AbstractParticles, nm::Symbol)
    val = p[nm]
    T = _get_eltype(val)
    if T <: Vector
        return eachcol(val)
    else
        return val isa AbstractVector ? val : fill(val, nrow(p))
    end
end

function Tables.getcolumn(p::AbstractParticles, i::Int)
    j = 0
    for key in Tables.columnnames(p)
        j += 1
        i == j && return Tables.getcolumn(p, key)
    end

    error("Column $(i) is out of bounds, only $j columns available.")
end

Tables.columnnames(p::AbstractParticles) = ParticlesColumnKeySet(p)

Tables.DataAPI.ncol(p::AbstractParticles) = length(Tables.columnnames(p))
Tables.DataAPI.nrow(p::AbstractParticles) = particle_number(p)


# iterate columns as vector accessors
struct ParticlesColumnKeySet{P<:AbstractParticles} <: AbstractSet{Symbol}
    particles::P
end

function Base.show(io::IO, mime::MIME"text/plain", pk::ParticlesColumnKeySet)
    print(io, "ParticlesColumnKeySet for a ")
    summary(io, pk.particles)
end

function Base.length(pk::ParticlesColumnKeySet)
    i = 0
    for _ in pk
        i += 1
    end
    return i
end

Base.isempty(pk::ParticlesColumnKeySet) = isempty(pk.particles)

function Base.iterate(pk::ParticlesColumnKeySet)
    isempty(pk) && return nothing

    p = pk.particles
    dkey, dkey_state = iterate(keys(p))

    vals = p[dkey]
    if vals isa AbstractMatrix
        dims = size(vals, 1)
        key_col = _get_current_matrix_column_key(vals, dkey, dims, dims)
        return key_col, (dkey, dkey_state, dims)
    end

    return dkey, (dkey, dkey_state, 1) # 1 is how many dimensions the dkey has
end

function Base.iterate(pk::ParticlesColumnKeySet, (dkey, dkey_state, state))
    p = pk.particles

    if state > 1
        vals = p[dkey]
        dims = size(vals, 1)
        state -= 1
        key_col = _get_current_matrix_column_key(vals, dkey, dims, state)
        return key_col, (dkey, dkey_state, state)
    end


    it = iterate(keys(p), dkey_state)
    isnothing(it) && return nothing

    dkey, dkey_state = it

    vals = p[dkey]
    if vals isa AbstractMatrix
        dims = size(vals, 1)
        key_col = _get_current_matrix_column_key(vals, dkey, dims, dims)
        return key_col, (dkey, dkey_state, dims)
    end

    return dkey, (dkey, dkey_state, 1) # 1 is how many dimensions the dkey has
end

function _get_current_matrix_column_key(vals::AbstractMatrix, dkey, dims, i)
    ind = dims - i + 1

    if vals isa NamedRowMatrix
        postfix = string("_", names(vals)[ind])
    else
        postfix = dims ≤ 3 ? ("_x", "_y", "_z")[ind] : "_$ind"
    end

    return Symbol(dkey, postfix)
end
