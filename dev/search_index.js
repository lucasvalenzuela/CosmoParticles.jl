var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = CosmoParticles","category":"page"},{"location":"#CosmoParticles","page":"Home","title":"CosmoParticles","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for CosmoParticles.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [CosmoParticles]","category":"page"},{"location":"#CosmoParticles.AbstractParticles","page":"Home","title":"CosmoParticles.AbstractParticles","text":"AbstractParticles\n\nAbstract supertype for storing particle data efficiently in arrays.\n\nParticles generally have IDs, positions, and velocities as properties. The properties are expected to be stored in a Dict{Symbol,Any} as vectors of length N, matrices of size mN, or as scalars (when the property is equal for all particles).\n\nThe inbuilt functionality of AbstractParticles includes accessing and setting properties via the following syntax:\n\np.id = [1, 2, 3]\np[:id] === p.id\n\nIf the struct has additional fields, these can also be accessed by p.field, but not by p[:field]. The latter syntax can only be used to access the particle properties.\n\nProperty keys\n\n:id: Vector of IDs\n:pos: 2N or 3N Matrix with positions\n:vel: 2N or 3N Matrix with velocities\n\nMethods\n\nThe methods Base.keys, Base.values, and Base.haskey are forwarded to the property Dict.\n\nConcrete types of AbstractParticles should have the following methods implemented (also see the implementation of Particles):\n\nBase.copy: returns new object containing a copy of the Dict\nBase.:(==)\nCosmoParticles.particle_name: returns the name of the struct to be printed via Base.show\nBase.propertynames: implement this if there are additional struct fields\nBase.show(io, mime, p)\n\n\n\n\n\n","category":"type"},{"location":"#CosmoParticles.Particles","page":"Home","title":"CosmoParticles.Particles","text":"struct Particles <: AbstractParticles\n    type::Symbol\n    props::Dict{Symbol,Any}\nend\n\nParticles of a certain type (typically something like :dm or :gas in the cosmological context) with their properties.\n\n\n\n\n\n","category":"type"},{"location":"#CosmoParticles.Particles-Tuple{Any}","page":"Home","title":"CosmoParticles.Particles","text":"Particles(type[, pairs::Pair...])\n\nCreate a Particles object with the given type and pairs of Symbol to the property values. These are passed to the underlying Dict. If no pairs are passed to the method, an empty Dict is created.\n\n\n\n\n\n","category":"method"},{"location":"#Base.filter!-Tuple{AbstractParticles}","page":"Home","title":"Base.filter!","text":"Base.filter!(p::AbstractParticles, ids)\n\nFilter the particles in-place by keeping only those with the given IDs.\n\n\n\n\n\n","category":"method"},{"location":"#Base.filter!-Tuple{Any, AbstractParticles}","page":"Home","title":"Base.filter!","text":"Base.filter!(f, p::AbstractParticles)\n\nFilter the particles in-place by a mask returned by the function f.\n\nThe function f takes the particles as argument and has to either return a BitArray mask with a length of the number of particles or an array of indices.\n\n\n\n\n\n","category":"method"},{"location":"#Base.filter-Tuple{AbstractParticles}","page":"Home","title":"Base.filter","text":"Base.filter(p::AbstractParticles, ids; affect=())\n\nCreate new particles with them filtered by keeping only those with the given IDs.\n\nIf the keyword argument affect is a non-empty tuple of Symbols, only those properties are filtered and added to the newly created particles object.\n\n\n\n\n\n","category":"method"},{"location":"#Base.filter-Tuple{Any, AbstractParticles}","page":"Home","title":"Base.filter","text":"Base.filter(f, p::AbstractParticles; affect=())\n\nCreate new particles with them filtered by a mask returned by the function f.\n\nThe function f takes the particles as argument and has to either return a BitArray mask with a length of the number of particles or an array of indices. If the keyword argument affect is a non-empty tuple of Symbols, only those properties are filtered and added to the newly created particles object.\n\n\n\n\n\n","category":"method"},{"location":"#Base.sort!-Tuple{AbstractParticles, Symbol}","page":"Home","title":"Base.sort!","text":"Base.sort!(p::AbstractParticles, prop::Symbol; kwargs...)\n\nSort the particles in-place according to the specified property.\n\nAll property arrays of the particles are rearranged in-place according to the property being sorted. The keyword arguments are passed on to Base.sortperm.\n\nThe sorting algorithm for unitful properties may also be RadixSort from SortingAlgorithms.jl.\n\n\n\n\n\n","category":"method"},{"location":"#Base.sort-Tuple{AbstractParticles, Symbol}","page":"Home","title":"Base.sort","text":"Base.sort(p::AbstractParticles, prop::Symbol; affect=(), kwargs...)\n\nCreate new particles with the particles sorted according to the specified property.\n\nAll property arrays of the particles are rearranged according to the property being sorted. If the keyword argument affect is a non-empty tuple of Symbols, only those properties are rearranged and added to the newly created particles object (besides the originally sorted property). The other keyword arguments are passed on to Base.sortperm.\n\nThe sorting algorithm for unitful properties may also be RadixSort from SortingAlgorithms.jl.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoParticles.get_props-Tuple{AbstractParticles}","page":"Home","title":"CosmoParticles.get_props","text":"CosmoParticles.get_props(p::AbstractParticles)\n\nReturn the property Dict belonging to the particles.\n\nThis returns p.props by default if not overridden.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoParticles.matrix_rotate!-Tuple{AbstractMatrix{<:Number}, AbstractMatrix{<:Real}}","page":"Home","title":"CosmoParticles.matrix_rotate!","text":"CosmoParticles.matrix_rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})\n\nIn-place version of CosmoParticles.matrix_rotate.\n\nThis is not exported.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoParticles.matrix_rotate-Tuple{AbstractMatrix{<:Number}, AbstractMatrix{<:Real}}","page":"Home","title":"CosmoParticles.matrix_rotate","text":"CosmoParticles.matrix_rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})\n\nReturns vals rotated by the rotation matrix rotmat by regular matrix multiplication. Works for any dimensions and is optimized for unitful arrays.\n\nThis is not exported.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoParticles.particle_name-Tuple{AbstractParticles}","page":"Home","title":"CosmoParticles.particle_name","text":"CosmoParticles.particle_name(::AbstractParticles)\n\nReturns the name of the particle type, which is used when printing an object of this type via Base.show.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoParticles.particle_number-Tuple{AbstractParticles}","page":"Home","title":"CosmoParticles.particle_number","text":"CosmoParticles.particle_number(p::AbstractParticles)\n\nReturns the number of particles by determining the length of one of the property arrays.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoParticles.rotate","page":"Home","title":"CosmoParticles.rotate","text":"rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop)\nrotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props=(:pos, :vel))\n\nRotates the specified properties props of the particles p by the rotation matrix rotmat.\n\nCreates a copy of the particles with only new pointers to the rotated properties (by default position and velocity). This means that mutations to other properties will affect both the original particles as well as the rotated copied particles. The properties should be given as a single Symbol, or an array or tuple of Symbols. Only the existing quantities will be rotated (e.g., this method does not throw an error if velocities are not given for the particles).\n\n\n\n\n\n","category":"function"},{"location":"#CosmoParticles.show_properties-Tuple{IO, AbstractString, AbstractParticles}","page":"Home","title":"CosmoParticles.show_properties","text":"CosmoParticles.show_properties(io::IO, mime, p::AbstractParticles)\n\nPrints the number of particles, the name of the type, and the property names.\n\nShould be used internally when overriding Base.show for concrete implementations of AbstractParticles.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoParticles.translate","page":"Home","title":"CosmoParticles.translate","text":"translate(p::AbstractParticles, Δx, prop::Symbol=:pos)\n\nTranslates the specified property of the particles p by Δx.\n\nCreates a copy of the particles with only new pointers to the translated property. Typically, Δx will be an AbstractVector, e.g., for positions or velocities.\n\n\n\n\n\n","category":"function"},{"location":"#CosmoParticles.translate!","page":"Home","title":"CosmoParticles.translate!","text":"translate!(p::AbstractParticles, Δx, prop::Symbol=:pos)\n\nIn-place version of translate.\n\n\n\n\n\n","category":"function"},{"location":"#LinearAlgebra.rotate!","page":"Home","title":"LinearAlgebra.rotate!","text":"LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props=(:pos, :vel))\nLinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop)\n\nIn-place version of rotate.\n\nThis function is reexported.\n\n\n\n\n\n","category":"function"}]
}
