Tables.istable(::Type{<:AbstractParticles}) = true
Tables.columnaccess(::Type{<:AbstractParticles}) = true
Tables.columns(p::AbstractParticles) = p

# The table has one column per vector property, per scalar property (repeated for every particle) and per row of a
# matrix property. The rows of a matrix are named as in the underscore access of `getindex`: by their row names for a
# `NamedRowMatrix`, by `x`, `y`, `z` for two or three rows, and by their numbers otherwise.
function _table_columns(p::AbstractParticles)
    cols = Tuple{Symbol,Symbol,Int}[] # column name, property, and row of a matrix property (0 for others)
    for key in keys(p)
        vals = p[key]
        if vals isa AbstractMatrix
            for (row, suffix) in zip(1:size(vals, 1), _row_suffixes(vals))
                push!(cols, (Symbol(key, "_", suffix), key, row))
            end
        else
            push!(cols, (key, key, 0))
        end
    end
    return cols
end

_row_suffixes(vals::NamedRowMatrix) = rownames(vals)
_row_suffixes(vals::AbstractMatrix) = 2 ≤ size(vals, 1) ≤ 3 ? (:x, :y, :z) : 1:size(vals, 1)

function _table_column(p::AbstractParticles, key::Symbol, row::Int)
    vals = p[key]
    row == 0 || return @view vals[row, :]
    return vals isa AbstractVector ? vals : Fill(vals, particle_number(p))
end

_table_column_eltype(vals::AbstractArray) = eltype(vals)
_table_column_eltype(val) = typeof(val)

Tables.columnnames(p::AbstractParticles) = [name for (name, _, _) in _table_columns(p)]

function Tables.schema(p::AbstractParticles)
    cols = _table_columns(p)
    names = [name for (name, _, _) in cols]
    types = [_table_column_eltype(p[key]) for (_, key, _) in cols]
    return Tables.Schema(names, types)
end

function Tables.getcolumn(p::AbstractParticles, nm::Symbol)
    cols = _table_columns(p)
    i = findfirst(col -> col[1] === nm, cols)
    if isnothing(i)
        matrix_cols = [name for (name, key, row) in cols if key === nm && row > 0]
        isempty(matrix_cols) || throw(
            ArgumentError(
                "the matrix property :$nm is not a table column, its rows are the columns " *
                join(repr.(matrix_cols), ", "),
            ),
        )
        throw(ArgumentError("there is no table column :$nm"))
    end
    _, key, row = cols[i]
    return _table_column(p, key, row)
end

function Tables.getcolumn(p::AbstractParticles, i::Int)
    _, key, row = _table_columns(p)[i]
    return _table_column(p, key, row)
end

Tables.getcolumn(p::AbstractParticles, ::Type, i::Int, ::Symbol) = Tables.getcolumn(p, i)

Tables.DataAPI.ncol(p::AbstractParticles) = length(_table_columns(p))
Tables.DataAPI.nrow(p::AbstractParticles) = particle_number(p)
