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
    dkey, dkey_state = iterate(keys(pk.particles))
    return _first_column_key(pk.particles, dkey, dkey_state)
end

function Base.iterate(pk::ParticlesColumnKeySet, (dkey, dkey_state, dims, state, names))
    if state > 1
        state -= 1
        return _matrix_column_key(dkey, names, dims, state), (dkey, dkey_state, dims, state, names)
    end

    it = iterate(keys(pk.particles), dkey_state)
    isnothing(it) && return nothing
    return _first_column_key(pk.particles, it...)
end

# The iteration state holds the number of rows of a matrix property, how many of its columns are left, and its row
# names (`nothing` if it has none), so the row names are only retrieved once per matrix.
function _first_column_key(p, dkey::Symbol, dkey_state)
    vals = p[dkey]
    vals isa AbstractMatrix || return dkey, (dkey, dkey_state, 1, 1, nothing)

    dims = size(vals, 1)
    names = vals isa NamedRowMatrix ? rownames(vals) : nothing
    return _matrix_column_key(dkey, names, dims, dims), (dkey, dkey_state, dims, dims, names)
end

function _matrix_column_key(dkey, names, dims, i)
    ind = dims - i + 1

    if isnothing(names)
        postfix = dims ≤ 3 ? ("_x", "_y", "_z")[ind] : "_$ind"
    else
        postfix = string("_", names[ind])
    end

    return Symbol(dkey, postfix)
end
