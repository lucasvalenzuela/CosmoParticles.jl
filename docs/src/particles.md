```@meta
CurrentModule = CosmoParticles
```

# Particles

A set of objects with different properties are represented by subtypes of `AbstractParticles`.
While not absolutely necessary, these kind of objects generally have IDs, positions, and velocity as properties.
In the context of cosmological simulations, this applies to particle and galaxy data.

For an example implementation of a concrete subtype of `AbstractParticles`, see the source code of [`Particles`](@ref).

```@docs
AbstractParticles
Particles
Particles(type::Symbol)
Particles(p::AllParticles, props=nothing)
CosmoParticles.get_props
CosmoParticles.particle_name
CosmoParticles.particle_number
CosmoParticles.show_properties(io::IO, mime::AbstractString, p::AbstractParticles)
```


## Tables

Particles implement the [Tables.jl](https://github.com/JuliaData/Tables.jl) interface, so they can be passed to any
package that reads tables, for example `DataFrame(p)` from [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) or
`CSV.write("dm.csv", p)` from [CSV.jl](https://github.com/JuliaData/CSV.jl).
Each particle is a row of the table. The columns are:
- one column for each vector property;
- one column for each row of a matrix property, named like the underscore access of the property: `pos_x`,
  `pos_y` and `pos_z` for two or three rows, by the row names for a
  [`NamedRowMatrix`](@extref NamedRowArrays.NamedRowMatrix), and `pos_1`, `pos_2`, … otherwise.

```jldoctest
julia> using CosmoParticles, FillArrays, Tables

julia> p = Particles(:dm, :id => [1, 2, 3], :pos => [0.0 1.0 2.0; 0.5 0.5 0.5; 1.0 0.0 -1.0], :mass => Fill(2.0, 3));

julia> Tables.columntable(p)
(id = [1, 2, 3], pos_x = [0.0, 1.0, 2.0], pos_y = [0.5, 0.5, 0.5], pos_z = [1.0, 0.0, -1.0], mass = Fill(2.0, 3))

julia> nrow(p), ncol(p)
(3, 5)
```

A table is read back into particles with [`Particles(type, table)`](@ref Particles(::Symbol, ::Any)), for example
from `CSV.File("dm.csv")`. The columns of the matrix rows are combined into matrices again:

```jldoctest
julia> using CosmoParticles, FillArrays, Tables

julia> p = Particles(:dm, :id => [1, 2, 3], :pos => [0.0 1.0 2.0; 0.5 0.5 0.5; 1.0 0.0 -1.0], :mass => Fill(2.0, 3));

julia> Particles(:dm, Tables.columntable(p)) == p
true
```
