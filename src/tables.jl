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

Tables.getcolumn(p::AbstractParticles, i::Int) = Tables.getcolumn(p, Tables.columnnames(p)[i])
Tables.columnnames(p::AbstractParticles) = collect(keys(p))

Tables.DataAPI.ncol(p::AbstractParticles) = length(Tables.columnnames(p))
Tables.DataAPI.nrow(p::AbstractParticles) = particle_number(p)
