Tables.istable(::Type{<:AbstractParticles}) = true
Tables.columnaccess(::Type{<:AbstractParticles}) = true

# The columns are determined once for a table conversion, so that consumers that read the columns one by one do
# not determine them again for each column.
function Tables.columns(p::AbstractParticles)
    cols = _table_columns(p)
    values = OrderedDict{Symbol,AbstractVector}(name => _table_column(p, key, row) for (name, key, row) in cols)
    return Tables.DictColumnTable(_table_schema(p, cols), values)
end

# The table has one column per vector property and per row of a matrix property. The rows of a matrix are named as
# in the underscore access of `getindex`: by their row names for a `NamedRowMatrix`, by `x`, `y`, `z` for two or
# three rows, and by their numbers otherwise. The column names must be unique, which fails, e.g., for a vector
# property `:pos_x` next to the matrix property `:pos`, whose first row is also named `pos_x`.
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
    _check_unique_columns(cols)
    return cols
end

function _check_unique_columns(cols)
    allunique(name for (name, _, _) in cols) && return nothing
    sources = Dict{Symbol,Tuple{Symbol,Int}}()
    for (name, key, row) in cols
        haskey(sources, name) && throw(
            ArgumentError(
                "the table column :$name is given by both $(_column_source(sources[name]...)) and " *
                _column_source(key, row),
            ),
        )
        sources[name] = (key, row)
    end
end

_column_source(key, row) = row == 0 ? "the property :$key" : "row $row of the matrix property :$key"

_row_suffixes(vals::NamedRowMatrix) = rownames(vals)
_row_suffixes(vals::AbstractMatrix) = 2 ≤ size(vals, 1) ≤ 3 ? (:x, :y, :z) : 1:size(vals, 1)

function _table_column(p::AbstractParticles, key::Symbol, row::Int)
    vals = p[key]
    row == 0 || return @view vals[row, :]
    return vals
end


Tables.columnnames(p::AbstractParticles) = [name for (name, _, _) in _table_columns(p)]

Tables.schema(p::AbstractParticles) = _table_schema(p, _table_columns(p))

function _table_schema(p::AbstractParticles, cols)
    names = [name for (name, _, _) in cols]
    types = [eltype(p[key]) for (_, key, _) in cols]
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

"""
    Particles(type::Symbol, table)

Create particles of the given `type` from a table that implements the
[Tables.jl](https://github.com/JuliaData/Tables.jl) interface.

Consecutive columns that share the name before their last underscore become the rows of a matrix property if their
suffixes are `x`, `y`, and possibly `z`, or `1`, `2`, ..., in this order. For example, the columns `pos_x`, `pos_y`,
and `pos_z` become the ``3×N`` matrix `pos`, but the columns `pos_y` and `pos_x` stay vectors. Only the suffix after the last underscore is considered, so the columns `pos_x_1` and
`pos_x_2` become the ``2×N`` matrix `pos_x`. Columns are not combined if another column has the name of the matrix.
All other columns become vector properties, which includes the rows of a
[`NamedRowMatrix`](@extref NamedRowArrays.NamedRowMatrix) in the table of particles.

# Examples
```jldoctest
julia> using Tables

julia> p = Particles(:dm, :id => [1, 2], :pos => [1.0 2.0; 3.0 4.0; 5.0 6.0]);

julia> t = Tables.columntable(p);

julia> keys(t)
(:id, :pos_x, :pos_y, :pos_z)

julia> Particles(:dm, t) == p
true
```
"""
function Particles(type::Symbol, table)
    Tables.istable(table) || throw(ArgumentError("a value of type $(typeof(table)) is not a table"))
    cols = Tables.columns(table)
    props = OrderedDict{Symbol,Any}()
    for (key, colnames) in _table_props(collect(Symbol, Tables.columnnames(cols)))
        if length(colnames) == 1 && colnames[1] === key
            props[key] = Tables.getcolumn(cols, key)
        else
            props[key] = stack((Tables.getcolumn(cols, name) for name in colnames); dims=1)
        end
    end
    return Particles(type, props)
end

# Pairs of the property name and its column names. Consecutive columns that share the name before their last underscore
# form the rows of a matrix if their suffixes are `x`, `y`, and possibly `z`, or `1`, `2`, ..., in this order.
function _table_props(names::Vector{Symbol})
    props = Pair{Symbol,Vector{Symbol}}[]
    i = firstindex(names)
    while i ≤ lastindex(names)
        key, suffix = _split_suffix(names[i])
        n = 1 # number of matrix rows
        if !isnothing(suffix) && !(key in names)
            while i + n ≤ lastindex(names)
                expected = _row_suffix(suffix, n + 1)
                (!isnothing(expected) && _split_suffix(names[i + n]) == (key, expected)) || break
                n += 1
            end
        end
        if n ≥ 2
            push!(props, key => names[i:(i + n - 1)])
        else
            push!(props, names[i] => [names[i]])
        end
        i += n
    end
    return props
end

# The name before the last underscore and the suffix after it, or `nothing`s if there is no such underscore.
function _split_suffix(name::Symbol)
    str = String(name)
    i = findlast('_', str)
    (isnothing(i) || i == firstindex(str) || i == lastindex(str)) && return nothing, nothing
    return Symbol(str[begin:prevind(str, i)]), str[nextind(str, i):end]
end

# The suffix of row `k` of a matrix whose first row has the suffix `first`, or `nothing` if there is no such row.
function _row_suffix(first::AbstractString, k::Integer)
    first == "x" && return k ≤ 3 ? ("x", "y", "z")[k] : nothing
    first == "1" && return string(k)
    return nothing
end
