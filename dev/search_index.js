var documenterSearchIndex = {"docs":
[{"location":"particles/","page":"Particles","title":"Particles","text":"CurrentModule = CosmoParticles","category":"page"},{"location":"particles/#Particles","page":"Particles","title":"Particles","text":"","category":"section"},{"location":"particles/","page":"Particles","title":"Particles","text":"A set of objects with different properties are represented by subtypes of AbstractParticles. While not absolutely necessary, these kind of objects generally have IDs, positions, and velocity as properties. In the context of cosmological simulations, this applies to particle and galaxy data.","category":"page"},{"location":"particles/","page":"Particles","title":"Particles","text":"For an example implementation of a concrete subtype of AbstractParticles, see the source code of Particles.","category":"page"},{"location":"particles/","page":"Particles","title":"Particles","text":"AbstractParticles\nParticles\nParticles(type[, pairs::Pair...])\nCosmoParticles.get_props\nCosmoParticles.particle_name\nCosmoParticles.particle_number\nCosmoParticles.show_properties","category":"page"},{"location":"particles/#CosmoParticles.AbstractParticles","page":"Particles","title":"CosmoParticles.AbstractParticles","text":"abstract type AbstractParticles end\n\nAbstract supertype for storing particle data efficiently in arrays.\n\nParticles generally have IDs, positions, and velocities as properties. The properties are expected to be stored in a Dict{Symbol,Any} as vectors of length N, matrices of size mN, or as scalars (when the property is equal for all particles).\n\nThe inbuilt functionality of AbstractParticles includes accessing and setting properties via the following syntax:\n\np.id = [1, 2, 3]\np[:id] === p.id\n\nIf the struct has additional fields, these can also be accessed by p.field, but not by p[:field]. The latter syntax can only be used to access the particle properties.\n\nProperty keys\n\n:id: vector of IDs\n:pos: 2N or 3N matrix with positions\n:vel: 2N or 3N matrix with velocities\n\nAny arrays may be of the type AbstractArray, provided the arrays are 1-indexed.\n\nMethods\n\nThe methods Base.keys, Base.values, Base.haskey, Base.empty, Base.empty!, and Base.isempty are forwarded to the property Dict.\n\nConcrete types of AbstractParticles should have the following methods implemented (also see the implementation of Particles):\n\nBase.copy: returns new object containing a copy of the Dict\nBase.empty: returns an identical object, but with an empty Dict\nBase.:(==)\nCosmoParticles.particle_name: returns the name of the struct to be printed via Base.show\nBase.propertynames: implement this if there are additional struct fields\nBase.show(io, mime, p)\n\n\n\n\n\n","category":"type"},{"location":"particles/#CosmoParticles.Particles","page":"Particles","title":"CosmoParticles.Particles","text":"struct Particles <: AbstractParticles\n    type::Symbol\n    props::Dict{Symbol,Any}\nend\n\nParticles of a certain type (typically something like :dm or :gas in the cosmological context) with their properties.\n\nThe following naming convention is to be used for specific properties:\n\n:id: vector of IDs\n:pos: 2N or 3N matrix with positions\n:vel: 2N or 3N matrix with velocities\n:mass: vector of masses or a scalar if all particles have the same mass\n:temp: vector of temperatures\n\n\n\n\n\n","category":"type"},{"location":"geometry/","page":"Geometry","title":"Geometry","text":"CurrentModule = CosmoParticles","category":"page"},{"location":"geometry/#Geometry","page":"Geometry","title":"Geometry","text":"","category":"section"},{"location":"geometry/","page":"Geometry","title":"Geometry","text":"To filter particles in a particular volume, such as a cube, a sphere, or a cylinder, geometries can be used. This package provides a variety of different geometries for multiple dimensions, including the following:","category":"page"},{"location":"geometry/","page":"Geometry","title":"Geometry","text":"Hyperrectangle\n3D Rectangular Cuboid\n2D Rectangle\nHypersphere\n3D Sphere\n2D Circle\nCylinder (standing cylinder aligned with the z axis and arbitrary orientation)","category":"page"},{"location":"geometry/","page":"Geometry","title":"Geometry","text":"AbstractCosmoGeometry\nCosmoParticles.geometry_enclosing_corners\nCosmoParticles.geometry_enclosing_center\nCosmoParticles.mask_in\nCosmoHyperrectangle\nCosmoCuboid\nCosmoCube\nCosmoRectangle\nCosmoSquare\nCosmoHypersphere\nCosmoSphere\nCosmoCircle\nCosmoCylinder\nCosmoStandingCylinder","category":"page"},{"location":"geometry/#CosmoParticles.AbstractCosmoGeometry","page":"Geometry","title":"CosmoParticles.AbstractCosmoGeometry","text":"abstract type AbstractCosmoGeometry end\n\nAbstract type for (multi-dimensional) geometry volumes, particularly for filtering with filter.\n\nAny subtypes of AbstractCosmoGeometry have to implement the following methods:\n\nCosmoParticles.geometry_enclosing_corners\nCosmoParticles.geometry_enclosing_center\nCosmoParticles.mask_in\n\n\n\n\n\n","category":"type"},{"location":"geometry/#CosmoParticles.geometry_enclosing_corners","page":"Geometry","title":"CosmoParticles.geometry_enclosing_corners","text":"geometry_enclosing_corners(geo::AbstractCosmoGeometry)\n\nReturn the lower left and upper right corners of the enclosing box of the geometry as a tuple of vectors.\n\nThe enclosing box is not necessarily the tightest fitting box.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"geometry/#CosmoParticles.geometry_enclosing_center","page":"Geometry","title":"CosmoParticles.geometry_enclosing_center","text":"geometry_enclosing_center(geo::AbstractCosmoGeometry)\n\nReturn the center of the enclosing box of the geometry as a vector.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"geometry/#CosmoParticles.mask_in","page":"Geometry","title":"CosmoParticles.mask_in","text":"mask_in(pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)\n\nReturn the BitArray mask of the positions (mathrmdims times N) located within the geometry.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"geometry/#CosmoParticles.CosmoHyperrectangle","page":"Geometry","title":"CosmoParticles.CosmoHyperrectangle","text":"struct CosmoHyperrectangle{T,N} <: AbstractCosmoGeometry where {T<:Number}\n    lowerleft::Vector{T}\n    upperright::Vector{T}\nend\n\nHyperrectangle aligned with the coordinate system axes given by its lower left and upper right corners.\n\nThe dimensions of space are given by N.\n\n\n\n\n\n","category":"type"},{"location":"geometry/#CosmoParticles.CosmoCuboid","page":"Geometry","title":"CosmoParticles.CosmoCuboid","text":"CosmoCuboid{T} = CosmoHyperrectangle{T,3}\n\nAlias for a 3D CosmoHyperrectangle.\n\n\n\n\n\n","category":"type"},{"location":"geometry/#CosmoParticles.CosmoCube","page":"Geometry","title":"CosmoParticles.CosmoCube","text":"CosmoCube(center::AbstractVector{<:Number}, radius::Number)\n\nReturn a cubic CosmoCuboid centered around center, with equal sidelengths.\n\n\n\n\n\n","category":"function"},{"location":"geometry/#CosmoParticles.CosmoRectangle","page":"Geometry","title":"CosmoParticles.CosmoRectangle","text":"CosmoRectangle{T} = CosmoHyperrectangle{T,2}\n\nAlias for a 2D CosmoHyperrectangle.\n\n\n\n\n\n","category":"type"},{"location":"geometry/#CosmoParticles.CosmoSquare","page":"Geometry","title":"CosmoParticles.CosmoSquare","text":"CosmoSquare(center::AbstractVector{<:Number}, radius::Number)\n\nReturn a square CosmoRectangle centered around center, with equal sidelengths.\n\n\n\n\n\n","category":"function"},{"location":"geometry/#CosmoParticles.CosmoHypersphere","page":"Geometry","title":"CosmoParticles.CosmoHypersphere","text":"struct CosmoHypersphere{T,N} <: AbstractCosmoGeometry where {T<:Number}\n    center::Vector{T}\n    radius::T\nend\n\nHypersphere given by its center and radius.\n\nIf different types are passed to the constructor, the types will be promoted without throwing an error.\n\n\n\n\n\n","category":"type"},{"location":"geometry/#CosmoParticles.CosmoSphere","page":"Geometry","title":"CosmoParticles.CosmoSphere","text":"CosmoSphere{T} = CosmoHypersphere{T,3}\n\nAlias for a 2D CosmoHypersphere.\n\n\n\n\n\n","category":"type"},{"location":"geometry/#CosmoParticles.CosmoCircle","page":"Geometry","title":"CosmoParticles.CosmoCircle","text":"CosmoCircle{T} = CosmoHypersphere{T,2}\n\nAlias for a 2D CosmoHypersphere.\n\n\n\n\n\n","category":"type"},{"location":"geometry/#CosmoParticles.CosmoCylinder","page":"Geometry","title":"CosmoParticles.CosmoCylinder","text":"struct CosmoCylinder{T} <: AbstractCosmoGeometry where {T<:Number}\n    startpos::Vector{T}\n    endpos::Vector{T}\n    radius::T\nend\n\nCylinder given by its end points startpos and endpos and radius.\n\nIf different types are passed to the constructor, the types will be promoted without throwing an error.\n\n\n\n\n\n","category":"type"},{"location":"geometry/#CosmoParticles.CosmoStandingCylinder","page":"Geometry","title":"CosmoParticles.CosmoStandingCylinder","text":"struct CosmoStandingCylinder{T} <: AbstractCosmoGeometry where {T<:Number}\n    center::Vector{T}\n    height::T\n    radius::T\nend\n\nStanding cylinder given by its center, height, and radius.\n\nStanding means that the cylinder is oriented such that its axis is aligned with the z axis. If different types are passed to the constructor, the types will be promoted without throwing an error.\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = CosmoParticles","category":"page"},{"location":"#CosmoParticles","page":"Home","title":"CosmoParticles","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CosmoParticles provides structs and a clean interface for working with sets of particles and galaxies that have individual properties, especially made for dealing with the data extracted from cosmological simulations.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is currently only available via the local registry CosmoSimsRegistry, which can be loaded from the Julia REPL with the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\npkg\"registry add https://github.com/lucasvalenzuela/CosmoSimsRegistry\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"This only needs to be done once per Julia installation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The latest version of the package is available for Julia 1.7 and newer versions and can be installed with:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"CosmoParticles\")","category":"page"},{"location":"operations/","page":"Operations","title":"Operations","text":"CurrentModule = CosmoParticles","category":"page"},{"location":"operations/#Operations","page":"Operations","title":"Operations","text":"","category":"section"},{"location":"operations/#Transformations","page":"Operations","title":"Transformations","text":"","category":"section"},{"location":"operations/","page":"Operations","title":"Operations","text":"rotate\nLinearAlgebra.rotate!\ntranslate\ntranslate!","category":"page"},{"location":"operations/#CosmoParticles.rotate","page":"Operations","title":"CosmoParticles.rotate","text":"rotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop)\nrotate(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props=(:pos, :vel))\n\nRotates the specified properties props of the particles p by the rotation matrix rotmat.\n\nCreates a copy of the particles with only new pointers to the rotated properties (by default position and velocity). This means that mutations to other properties will affect both the original particles as well as the rotated copied particles. The properties should be given as a single Symbol, or an array or tuple of Symbols. Only the existing quantities will be rotated (e.g., this method does not throw an error if velocities are not given for the particles).\n\n\n\n\n\n","category":"function"},{"location":"operations/#LinearAlgebra.rotate!","page":"Operations","title":"LinearAlgebra.rotate!","text":"LinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, props=(:pos, :vel))\nLinearAlgebra.rotate!(p::AbstractParticles, rotmat::AbstractMatrix{<:Real}, prop)\n\nIn-place version of rotate.\n\nThis function is reexported.\n\n\n\n\n\n","category":"function"},{"location":"operations/#CosmoParticles.translate","page":"Operations","title":"CosmoParticles.translate","text":"translate(p::AbstractParticles, Δx, prop::Symbol=:pos)\n\nTranslates the specified property of the particles p by Δx.\n\nCreates a copy of the particles with only new pointers to the translated property. Typically, Δx will be an AbstractVector, e.g., for positions or velocities.\n\n\n\n\n\n","category":"function"},{"location":"operations/#CosmoParticles.translate!","page":"Operations","title":"CosmoParticles.translate!","text":"translate!(p::AbstractParticles, Δx, prop::Symbol=:pos)\n\nIn-place version of translate.\n\n\n\n\n\n","category":"function"},{"location":"operations/#Further-Operations","page":"Operations","title":"Further Operations","text":"","category":"section"},{"location":"operations/","page":"Operations","title":"Operations","text":"Base.sort\nBase.sort!\nBase.filter\nBase.filter!\nCosmoParticles.applyind\nCosmoParticles.applyind!\nCosmoParticles.findall_in\nCosmoParticles.findall_in_sorted","category":"page"},{"location":"operations/#Base.sort","page":"Operations","title":"Base.sort","text":"Base.sort(p::AbstractParticles, prop::Symbol; affect=keys(p), kwargs...)\n\nCreate new particles with the particles sorted according to the specified property.\n\nAll property arrays of the particles are rearranged according to the property being sorted. If the keyword argument affect is a non-empty tuple of Symbols, only those properties are rearranged and added to the newly created particles object. The other keyword arguments are passed on to Base.sortperm.\n\nThe sorting algorithm for unitful properties may also be RadixSort from SortingAlgorithms.jl.\n\n\n\n\n\n","category":"function"},{"location":"operations/#Base.sort!","page":"Operations","title":"Base.sort!","text":"Base.sort!(p::AbstractParticles, prop::Symbol; kwargs...)\n\nSort the particles in-place according to the specified property.\n\nAll property arrays of the particles are rearranged in-place according to the property being sorted. The keyword arguments are passed on to Base.sortperm.\n\nThe sorting algorithm for unitful properties may also be RadixSort from SortingAlgorithms.jl.\n\n\n\n\n\n","category":"function"},{"location":"operations/#Base.filter","page":"Operations","title":"Base.filter","text":"Base.filter(f, p::AbstractParticles; affect=keys(p))\n\nCreate new particles with them filtered by a mask returned by the function f.\n\nThe function f takes the particles as argument and has to either return a BitArray mask with a length of the number of particles or an array of indices. If the keyword argument affect is a non-empty tuple of Symbols, only those properties are filtered and added to the newly created particles object.\n\n\n\n\n\nBase.filter(p::AbstractParticles, ids; affect=keys(p))\n\nCreate new particles with them filtered by keeping only those with the given IDs.\n\nIf the keyword argument affect is a non-empty tuple of Symbols, only those properties are filtered and added to the newly created particles object.\n\n\n\n\n\nBase.filter(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos; affect=keys(p))\n\nCreate new particles with them filtered by keeping only those inside the given AbstractCosmoGeometry.\n\nThe filter is applied to the property specified. If the keyword argument affect is a non-empty tuple of Symbols, only those properties are filtered and added to the newly created particles object.\n\n\n\n\n\n","category":"function"},{"location":"operations/#Base.filter!","page":"Operations","title":"Base.filter!","text":"Base.filter!(f, p::AbstractParticles)\n\nFilter the particles in-place by a mask returned by the function f.\n\nThe function f takes the particles as argument and has to either return a BitArray mask with a length of the number of particles or an array of indices.\n\n\n\n\n\nBase.filter!(p::AbstractParticles, ids)\n\nFilter the particles in-place by keeping only those with the given IDs.\n\n\n\n\n\nBase.filter!(p::AbstractParticles, geo::AbstractCosmoGeometry, prop::Symbol=:pos)\n\nFilter the particles in-place by keeping only those inside the given AbstractCosmoGeometry.\n\nThe filter is applied to the property specified.\n\n\n\n\n\n","category":"function"},{"location":"operations/#CosmoParticles.applyind","page":"Operations","title":"CosmoParticles.applyind","text":"CosmoParticles.applyind(p::AbstractParticles, ind::AbstractVector; affect=keys(p))\n\nCreate new particles with the given indices or mask applied to all particle properties.\n\nThis can also be called by the simple syntax p[ind]. If the keyword argument affect is a non-empty tuple of Symbols, only those properties are indexed into and added to the newly created particles object.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"operations/#CosmoParticles.applyind!","page":"Operations","title":"CosmoParticles.applyind!","text":"CosmoParticles.applyind!(p::AbstractParticles, ind::AbstractVector)\n\nIn-place application of indices or mask to all particle properties.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"operations/#CosmoParticles.findall_in","page":"Operations","title":"CosmoParticles.findall_in","text":"CosmoParticles.findall_in(a::AbstractVector, set)\n\nReturn all indices of a that are in set.\n\nIf both a and set are sorted AbstractVectors, then the optimized findall_in_sorted is called. Otherwise, a Set is constructed from the Vector to perform the checks with in.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"operations/#CosmoParticles.findall_in_sorted","page":"Operations","title":"CosmoParticles.findall_in_sorted","text":"CosmoParticles.findall_in_sorted(a::AbstractVector, set::AbstractVector)\n\nReturn all indices of a that are in set, where both a and set are assumed to be sorted.\n\nThis uses an optimized algorithm that is faster than creating a Set from set and performing checks with in.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"operations/#Internals","page":"Operations","title":"Internals","text":"","category":"section"},{"location":"operations/","page":"Operations","title":"Operations","text":"CosmoParticles.matrix_rotate\nCosmoParticles.matrix_rotate!\nCosmoParticles._applyind","category":"page"},{"location":"operations/#CosmoParticles.matrix_rotate","page":"Operations","title":"CosmoParticles.matrix_rotate","text":"CosmoParticles.matrix_rotate(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})\n\nReturns vals rotated by the rotation matrix rotmat by regular matrix multiplication. Works for any dimensions and is optimized for unitful arrays.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"operations/#CosmoParticles.matrix_rotate!","page":"Operations","title":"CosmoParticles.matrix_rotate!","text":"CosmoParticles.matrix_rotate!(vals::AbstractMatrix{<:Number}, rotmat::AbstractMatrix{<:Real})\n\nIn-place version of CosmoParticles.matrix_rotate.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"},{"location":"operations/#CosmoParticles._applyind","page":"Operations","title":"CosmoParticles._applyind","text":"CosmoParticles._applyind(a, ind::AbstractVector)\n\nApply indices or mask to a Number, Vector, or Matrix.\n\nThe following indexing is applied:\n\na::Number: a is returned directly.\na::Vector: a[ind] is returned.\na::Matrix: a[:, ind] is returned.\n\nThis is not exported.\n\n\n\n\n\n","category":"function"}]
}
