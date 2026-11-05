"""
    abstract type AbstractCosmoGeometry end

Abstract type for (multi-dimensional) geometry volumes, particularly for filtering with [`filter`](@ref).

Any subtypes of `AbstractCosmoGeometry` have to implement the following methods:
- [`CosmoParticles.geometry_enclosing_corners`](@ref)
- [`CosmoParticles.mask_in!`](@ref)

The following methods are optional, but have default implementations:
- [`CosmoParticles.geometry_enclosing_center`](@ref)
- [`CosmoParticles.mask_in`](@ref)
- [`CosmoParticles.translate`](@ref)

# Examples

```jldoctest
julia> p = Particles(:dm, :id => [1, 2, 3], :pos => [0.5 2.0 0.0; 0.0 0.0 -0.5; 0.0 0.0 0.0]);

julia> filter(p, CosmoSphere(1.0)).id
2-element Vector{Int64}:
 1
 3
```
"""
abstract type AbstractCosmoGeometry end

"""
    CosmoParticles.geometry_enclosing_corners(geo::AbstractCosmoGeometry)

Return the lower left and upper right corners of the enclosing box of the geometry as a tuple of vectors.

The enclosing box is not necessarily the tightest fitting box.

# Examples

```jldoctest
julia> s = CosmoSphere([1.0, 2.0, 3.0], 1.0);

julia> CosmoParticles.geometry_enclosing_corners(s)
([0.0, 1.0, 2.0], [2.0, 3.0, 4.0])
```
"""
function geometry_enclosing_corners end

"""
    CosmoParticles.geometry_enclosing_center(geo::AbstractCosmoGeometry)

Return the center of the enclosing box of the geometry as a vector.

# Examples

```jldoctest
julia> c = CosmoCuboid([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);

julia> CosmoParticles.geometry_enclosing_center(c)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 0.5
 1.0
 1.5
```
"""
function geometry_enclosing_center(geo::AbstractCosmoGeometry)
    corners = geometry_enclosing_corners(geo)
    return 1 // 2 .* (corners[1] .+ corners[2])
end

@doc raw"""
    CosmoParticles.mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)

Return the `BitVector` mask of the positions (``\mathrm{dims} × N``) located within the geometry.

# Examples

```jldoctest
julia> s = CosmoSphere(1.0);

julia> pos = [0.5 2.0 5.0; 0.0 0.0 0.0; 0.0 0.0 0.0];

julia> CosmoParticles.mask_in!(BitVector(undef, 3), pos, s)
3-element BitVector:
 1
 0
 0
```
"""
function mask_in! end

@doc raw"""
    CosmoParticles.mask_in(pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)

Return the `BitArray` mask of the positions (``\mathrm{dims} × N``) located within the geometry.

Calls [`mask_in!`](@ref) by default.

# Examples

```jldoctest
julia> s = CosmoSphere(1.0);

julia> pos = [0.5 2.0 5.0; 0.0 0.0 0.0; 0.0 0.0 0.0];

julia> CosmoParticles.mask_in(pos, s)
3-element BitVector:
 1
 0
 0
```
"""
function mask_in(pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)
    mask = BitVector(undef, size(pos, 2))
    return mask_in!(mask, pos, geo)
end

function _check_mask_dims(mask::AbstractVector, pos::AbstractMatrix, N)
    size(pos, 1) == N || throw(
        DimensionMismatch("the geometry is $N-dimensional, but the positions have $(size(pos, 1)) rows"),
    )
    size(pos, 2) == length(mask) || throw(
        DimensionMismatch("the mask has length $(length(mask)), but there are $(size(pos, 2)) positions"),
    )
    return nothing
end

# The geometries combined by a union, intersection, or difference must have the same number of dimensions.
function _check_geometries(geos::Tuple)
    isempty(geos) && throw(ArgumentError("at least one geometry is required"))
    for geo in geos
        geo isa AbstractCosmoGeometry || throw(ArgumentError("expected geometries, got a $(typeof(geo))"))
    end
    dims = [length(first(geometry_enclosing_corners(geo))) for geo in geos]
    allequal(dims) || throw(
        DimensionMismatch("the geometries must have the same number of dimensions, got $(join(dims, ", "))"),
    )
    return nothing
end

# Adding zero turns -0.0 into 0.0, which `==` treats as equal to 0.0 but `hash` does not.
_hash_normalized(x::Union{Number,AbstractArray{<:Number}}, h::UInt) = hash(x .+ zero.(x), h)
_hash_normalized(x, h::UInt) = hash(x, h)

# For geometries that are determined by their fields: equal if all fields are equal, with a matching hash.
_fields_equal(g1, g2) = all(f -> getfield(g1, f) == getfield(g2, f), fieldnames(typeof(g1)))
function _fields_hash(g, h::UInt)
    return foldl(
        (h, f) -> _hash_normalized(getfield(g, f), h),
        fieldnames(typeof(g));
        init=hash(nameof(typeof(g)), h),
    )
end



"""
    struct CosmoUnionGeometry <: AbstractCosmoGeometry
        geos::Tuple
    end

Union of multiple geometries.

This union geometry represents all of the volumes contained within the geometries. 

This type is not exported; it is accessible as `CosmoParticles.CosmoUnionGeometry`.
Union geometries are created with [`union`](@ref Base.union(::AbstractCosmoGeometry...)).
"""
struct CosmoUnionGeometry <: AbstractCosmoGeometry
    geos::Tuple

    function CosmoUnionGeometry(geos)
        geos = Tuple(geos)
        _check_geometries(geos)
        new(geos)
    end
end

"""
    CosmoParticles.CosmoUnionGeometry(geos::AbstractVector)
    CosmoParticles.CosmoUnionGeometry(geos::AbstractCosmoGeometry...)

Create a union of geometries from the passed geometries.
"""
CosmoUnionGeometry(geos::AbstractCosmoGeometry...) = CosmoUnionGeometry(geos)

"""
    Base.union(geos::AbstractCosmoGeometry...) 

Create a union of geometries from the passed geometries.

# Examples

```jldoctest
julia> g1 = CosmoSphere(1.0);

julia> g2 = CosmoSphere([3.0, 0.0, 0.0], 1.0);

julia> geo = union(g1, g2);

julia> pos = [0.5 2.0 5.0; 0.0 0.0 0.0; 0.0 0.0 0.0];

julia> CosmoParticles.mask_in(pos, geo) == CosmoParticles.mask_in(pos, g1) .| CosmoParticles.mask_in(pos, g2)
true
```
"""
Base.union(geos::AbstractCosmoGeometry...) = CosmoUnionGeometry(geos)

function geometry_enclosing_corners(g::CosmoUnionGeometry)
    corners = [geometry_enclosing_corners(geo) for geo in g.geos]
    lowerleft = reduce((x, y) -> min.(x, y), (ll for (ll, _) in corners))
    upperright = reduce((x, y) -> max.(x, y), (ur for (_, ur) in corners))

    return lowerleft, upperright
end

function mask_in(pos::AbstractMatrix{<:Number}, g::CosmoUnionGeometry)
    reduce(.|, mask_in(pos, geo) for geo in g.geos)
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, g::CosmoUnionGeometry)
    mask .= mask_in(pos, g)
end

function translate(g::CosmoUnionGeometry, Δx::AbstractVector{<:Number})
    CosmoUnionGeometry(Tuple(translate(geo, Δx) for geo in g.geos))
end


"""
    struct CosmoIntersectGeometry <: AbstractCosmoGeometry
        geos::Tuple
    end

Intersect of multiple geometries.

This intersect geometry represents the volume that is shared among all of the geometries. 

This type is not exported; it is accessible as `CosmoParticles.CosmoIntersectGeometry`.
Intersect geometries are created with [`intersect`](@ref Base.intersect(::AbstractCosmoGeometry...)).
"""
struct CosmoIntersectGeometry <: AbstractCosmoGeometry
    geos::Tuple

    function CosmoIntersectGeometry(geos)
        geos = Tuple(geos)
        _check_geometries(geos)
        new(geos)
    end
end

"""
    CosmoParticles.CosmoIntersectGeometry(geos::AbstractVector)
    CosmoParticles.CosmoIntersectGeometry(geos::AbstractCosmoGeometry...)

Create an intersect of geometries from the passed geometries.
"""
CosmoIntersectGeometry(geos::AbstractCosmoGeometry...) = CosmoIntersectGeometry(geos)

"""
    Base.intersect(geos::AbstractCosmoGeometry...) 

Create an intersect of geometries from the passed geometries.

# Examples

```jldoctest
julia> g1 = CosmoSphere(2.0);

julia> g2 = CosmoSphere([3.0, 0.0, 0.0], 2.0);

julia> geo = intersect(g1, g2);

julia> pos = [0.5 1.5 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0];

julia> CosmoParticles.mask_in(pos, geo) == CosmoParticles.mask_in(pos, g1) .& CosmoParticles.mask_in(pos, g2)
true
```
"""
Base.intersect(geos::AbstractCosmoGeometry...) = CosmoIntersectGeometry(geos)

function geometry_enclosing_corners(g::CosmoIntersectGeometry)
    corners = [geometry_enclosing_corners(geo) for geo in g.geos]
    lowerleft = reduce((x, y) -> max.(x, y), (ll for (ll, _) in corners))
    upperright = reduce((x, y) -> min.(x, y), (ur for (_, ur) in corners))

    if any(lowerleft .≥ upperright)
        return zero.(lowerleft), zero.(upperright)
    end

    return lowerleft, upperright
end

function mask_in(pos::AbstractMatrix{<:Number}, g::CosmoIntersectGeometry)
    reduce(.&, mask_in(pos, geo) for geo in g.geos)
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, g::CosmoIntersectGeometry)
    mask .= mask_in(pos, g)
end

function translate(g::CosmoIntersectGeometry, Δx::AbstractVector{<:Number})
    CosmoIntersectGeometry(Tuple(translate(geo, Δx) for geo in g.geos))
end


"""
    struct CosmoDiffGeometry <: AbstractCosmoGeometry
        geo::AbstractCosmoGeometry
        geos::Tuple
    end

Difference of a geometry with multiple geometries.

This difference geometry represents the volume of `geo`, but that is not contained in any of the geometries `geos`.

This type is not exported; it is accessible as `CosmoParticles.CosmoDiffGeometry`.
Difference geometries are created with [`setdiff`](@ref Base.setdiff(::AbstractCosmoGeometry, ::AbstractCosmoGeometry...)).
"""
struct CosmoDiffGeometry <: AbstractCosmoGeometry
    geo::AbstractCosmoGeometry
    geos::Tuple

    function CosmoDiffGeometry(geo::AbstractCosmoGeometry, geos)
        geos = Tuple(geos)
        isempty(geos) && throw(ArgumentError("at least one geometry to subtract is required"))
        _check_geometries((geo, geos...))
        new(geo, geos)
    end
end

"""
    CosmoParticles.CosmoDiffGeometry(geo::AbstractCosmoGeometry, geos::AbstractVector)
    CosmoParticles.CosmoDiffGeometry(geo::AbstractCosmoGeometry, geos::AbstractCosmoGeometry...)

Create a difference of geometries from the passed geometries.
"""
CosmoDiffGeometry(geo::AbstractCosmoGeometry, geos::AbstractCosmoGeometry...) = CosmoDiffGeometry(geo, geos)

"""
    Base.setdiff(geo::AbstractCosmoGeometry, geos::AbstractCosmoGeometry...) 

Create a difference of geometries from the passed geometries.

# Examples

```jldoctest
julia> g1 = CosmoSphere(2.0);

julia> g2 = CosmoSphere([1.0, 0.0, 0.0], 1.0);

julia> geo = setdiff(g1, g2);

julia> pos = [-1.0 0.5 1.5; 0.0 0.0 0.0; 0.0 0.0 0.0];

julia> CosmoParticles.mask_in(pos, geo) == CosmoParticles.mask_in(pos, g1) .& .~CosmoParticles.mask_in(pos, g2)
true
```
"""
Base.setdiff(geo::AbstractCosmoGeometry, geos::AbstractCosmoGeometry...) = CosmoDiffGeometry(geo, geos)

geometry_enclosing_corners(g::CosmoDiffGeometry) = geometry_enclosing_corners(g.geo)

function mask_in(pos::AbstractMatrix{<:Number}, g::CosmoDiffGeometry)
    mask = mask_in(pos, g.geo)
    masknot = reduce(.|, mask_in(pos, geo) for geo in g.geos)

    return @. mask & ~masknot
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, g::CosmoDiffGeometry)
    mask .= mask_in(pos, g)
end

function translate(g::CosmoDiffGeometry, Δx::AbstractVector{<:Number})
    CosmoDiffGeometry(translate(g.geo, Δx), Tuple(translate(geo, Δx) for geo in g.geos))
end


"""
    struct Rotated{G<:AbstractCosmoGeometry,R<:AbstractMatrix} <: AbstractCosmoGeometry
        geo::G
        rotmat::R
    end

Rotated geometry defined by a rotation matrix.

This can be created by [`rotate`](@ref) applied to a geometry.

This type is not exported; it is accessible as `CosmoParticles.Rotated`.

# Examples

```jldoctest
julia> R = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0];

julia> rot = rotate(CosmoCuboid([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0]), R);

julia> rot isa CosmoParticles.Rotated
true
```
"""
struct Rotated{G<:AbstractCosmoGeometry,R<:AbstractMatrix} <: AbstractCosmoGeometry
    geo::G
    rotmat::R
end

"""
    CosmoParticles.rotation_matrix(rot::Rotated)

Returns the rotation matrix that rotates the original geometry into the rotated state.

# Examples

```jldoctest
julia> R = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0];

julia> rot = rotate(CosmoCuboid([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0]), R);

julia> CosmoParticles.rotation_matrix(rot)
3×3 Matrix{Float64}:
 0.0  -1.0  0.0
 1.0   0.0  0.0
 0.0   0.0  1.0
```
"""
rotation_matrix(rot::Rotated) = rot.rotmat

"""
    CosmoParticles.rotation_matrix_inv(rot::Rotated)

Returns the rotation matrix that rotates the rotated geometry back to its original orientation.

Use this to rotate objects in the rotated frame of reference to the original frame of reference
of the geometry.

# Examples

```jldoctest
julia> R = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0];

julia> rot = rotate(CosmoCuboid([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0]), R);

julia> CosmoParticles.rotation_matrix_inv(rot) * CosmoParticles.rotation_matrix(rot)
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
```
"""
rotation_matrix_inv(rot::Rotated) = transpose(rot.rotmat)

function geometry_enclosing_corners(rot::Rotated)
    bounding_corners = geometry_enclosing_corners(rot.geo)
    n = length(bounding_corners[1])

    # create matrix of positions of the hypercuboid corners
    pos = Matrix{eltype(bounding_corners[1])}(undef, n, 2^n)
    @inbounds for bitval in axes(pos, 2)
        for i in 0:(n - 1)
            # bitval defines which corner vector it should take the corners from,
            # e.g. bitval = 1 for 3D space represents 100,
            # i.e. [upperright[1], lowerleft[2], lowerleft[3]]
            # or bitval = 6 for 3D space represents 110,
            # i.e. [upperright[1], upperright[2], lowerleft[3]]
            ind_corner = ifelse(iszero(((1 << i) & bitval)), 1, 2)
            pos[i + 1, bitval] = bounding_corners[ind_corner][i + 1]
        end
    end

    # rotate the hypercuboid corners by the matrix
    posrot = rot.rotmat * pos

    # extract the extrema along each dimension as the lower left and upper right corners
    lowerleft, upperright = vec(minimum(posrot; dims=2)), vec(maximum(posrot; dims=2))
    return _vector_like(bounding_corners[1], lowerleft), _vector_like(bounding_corners[2], upperright)
end

# The corners of a geometry are static vectors if those of the geometry it is derived from are.
_vector_like(::SVector{N}, v::AbstractVector) where {N} = SVector{N}(v)
_vector_like(_, v::AbstractVector) = v

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, rot::Rotated)
    mask_in!(mask, matrix_rotate(pos, rotation_matrix_inv(rot)), rot.geo)
end

"""
    rotate(geo::AbstractCosmoGeometry, rotmat::AbstractMatrix{<:Real})

Returns a rotated geometry rotated by the rotation matrix.

Individual geometries may specify that this returns a different type than [`Rotated`](@ref).

# Examples

```jldoctest
julia> R = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0];

julia> rot = rotate(CosmoCuboid([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0]), R);

julia> CosmoParticles.geometry_enclosing_corners(rot)
([-2.0, -1.0, -3.0], [2.0, 1.0, 3.0])
```
"""
rotate(geo::AbstractCosmoGeometry, rotmat::AbstractMatrix{<:Real}) = Rotated(geo, rotmat)
rotate(rot::Rotated, rotmat::AbstractMatrix{<:Real}) = Rotated(rot.geo, rotmat * rot.rotmat)



"""
    struct Translated{G<:AbstractCosmoGeometry,L<:AbstractVector} <: AbstractCosmoGeometry
        geo::G
        Δx::L
    end

Translated geometry defined by a translation vector.

This can be created by [`translate`](@ref) applied to a geometry.

This type is not exported; it is accessible as `CosmoParticles.Translated`.

# Examples

```jldoctest
julia> t = translate(CosmoEllipse(2.0; q=0.5), [1.0, 0.0]);

julia> t isa CosmoParticles.Translated
true

julia> CosmoParticles.geometry_enclosing_center(t)
2-element StaticArraysCore.SVector{2, Float64} with indices SOneTo(2):
 1.0
 0.0
```
"""
struct Translated{G<:AbstractCosmoGeometry,L<:AbstractVector} <: AbstractCosmoGeometry
    geo::G
    Δx::L
end

function geometry_enclosing_corners(tra::Translated)
    lowerleft, upperright = geometry_enclosing_corners(tra.geo)
    return lowerleft .+ tra.Δx, upperright .+ tra.Δx
end

function geometry_enclosing_center(tra::Translated)
    geometry_enclosing_center(tra.geo) .+ tra.Δx
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, tra::Translated)
    mask_in!(mask, pos .- tra.Δx, tra.geo)
end

"""
    translate(geo::AbstractCosmoGeometry, Δx::AbstractVector{<:Number})

Returns a new geometry translated by `Δx`.

Individual geometries may specify that this returns a different type than [`Translated`](@ref).

# Examples

```jldoctest
julia> t = translate(CosmoEllipse(2.0; q=0.5), [1.0, 0.0]);

julia> CosmoParticles.geometry_enclosing_corners(t)
([-1.0, -1.0], [3.0, 1.0])
```
"""
translate(geo::AbstractCosmoGeometry, Δx::AbstractVector{<:Number}) = Translated(geo, SVector{length(Δx)}(Δx))
translate(tra::Translated, Δx::AbstractVector{<:Number}) = Translated(tra.geo, tra.Δx .+ Δx)





"""
    struct CosmoHyperrectangle{T<:Number,N} <: AbstractCosmoGeometry
        lowerleft::SVector{N,T}
        upperright::SVector{N,T}
    end

`N`-dimensional hyperrectangle aligned with the coordinate system axes given by its lower left and upper right corners.

The dimensions of space are given by `N`.

# Examples

```jldoctest
julia> r = CosmoHyperrectangle([0.0, 0.0], [1.0, 2.0]);

julia> CosmoParticles.mask_in([0.5 2.0; 1.0 1.0], r)
2-element BitVector:
 1
 0
```
"""
struct CosmoHyperrectangle{T<:Number,N} <: AbstractCosmoGeometry
    lowerleft::SVector{N,T}
    upperright::SVector{N,T}

    function CosmoHyperrectangle{T,N}(
        lowerleft::AbstractVector{<:Number},
        upperright::AbstractVector{<:Number},
    ) where {T<:Number,N}
        length(lowerleft) == length(upperright) == N || throw(
            DimensionMismatch(
                "the corners must have length $N, got lengths $(length(lowerleft)) and $(length(upperright))",
            ),
        )
        new{T,N}(SVector{N,T}(lowerleft), SVector{N,T}(upperright))
    end
end

"""
    CosmoHyperrectangle(center::AbstractVector{<:Number}, sidelengths::Number...)
    CosmoHyperrectangle(center::AbstractVector{<:Number}; halflengths)

Return a `CosmoHyperrectangle` centered around `center`, with the given full `sidelengths` or `halflengths`.

The lengths are given in order of the axes, i.e., x, y or x, y, z for 2D and 3D, respectively.

# Examples

```jldoctest
julia> CosmoHyperrectangle([0.0, 0.0, 0.0, 0.0], 2, 4, 6, 8)
CosmoHyperrectangle([-1.0, -2.0, -3.0, -4.0], [1.0, 2.0, 3.0, 4.0])
```
"""
function CosmoHyperrectangle(center::AbstractVector{<:Number}, sidelength::Number, sidelengths::Number...)
    halflengths = [sidelength, sidelengths...]
    halflengths *= 1 // 2
    return CosmoHyperrectangle(center .- halflengths, center .+ halflengths)
end

function CosmoHyperrectangle(center::AbstractVector{<:Number}; halflengths)
    _check_halflengths(center, halflengths)
    return CosmoHyperrectangle(center .- halflengths, center .+ halflengths)
end

function _check_halflengths(center, halflengths)
    length(halflengths) == length(center) || throw(
        DimensionMismatch(
            "the half lengths must have length $(length(center)) like the center, got length $(length(halflengths))",
        ),
    )
    return nothing
end

function Base.:(==)(r1::CosmoHyperrectangle, r2::CosmoHyperrectangle)
    return r1.lowerleft == r2.lowerleft && r1.upperright == r2.upperright
end

function Base.hash(r::CosmoHyperrectangle, h::UInt)
    return _hash_normalized(r.upperright, _hash_normalized(r.lowerleft, hash(:CosmoHyperrectangle, h)))
end


"""
    CosmoCuboid{T} = CosmoHyperrectangle{T,3}

Alias for a 3D [`CosmoHyperrectangle`](@ref).
"""
const CosmoCuboid{T} = CosmoHyperrectangle{T,3} where {T<:Number}

"""
    CosmoCuboid(center::AbstractVector{<:Number}, x::Number, y::Number, z::Number)
    CosmoCuboid(center::AbstractVector{<:Number}; halflengths)

Return a `CosmoCuboid` centered around `center`, with the full side lengths `x`, `y`, and `z` or the three
`halflengths`.

# Examples

```jldoctest
julia> CosmoCuboid([0.0, 0.0, 0.0], 2, 4, 6)
CosmoCuboid([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0])
```
"""
function CosmoCuboid(center::AbstractVector{<:Number}, x::Number, y::Number, z::Number)
    halflengths = [x, y, z]
    halflengths *= 1 // 2
    return CosmoCuboid(center .- halflengths, center .+ halflengths)
end

function CosmoCuboid(center::AbstractVector{<:Number}; halflengths)
    _check_halflengths(center, halflengths)
    return CosmoCuboid(center .- halflengths, center .+ halflengths)
end


"""
    CosmoRectangle{T} = CosmoHyperrectangle{T,2}

Alias for a 2D [`CosmoHyperrectangle`](@ref).
"""
const CosmoRectangle{T} = CosmoHyperrectangle{T,2} where {T<:Number}

"""
    CosmoRectangle(center::AbstractVector{<:Number}, x::Number, y::Number)
    CosmoRectangle(center::AbstractVector{<:Number}; halflengths)

Return a `CosmoRectangle` centered around `center`, with the full side lengths `x` and `y` or the two
`halflengths`.

# Examples

```jldoctest
julia> CosmoRectangle([0.0, 0.0], 2, 4)
CosmoRectangle([-1.0, -2.0], [1.0, 2.0])
```
"""
function CosmoRectangle(center::AbstractVector{<:Number}, x::Number, y::Number)
    halflengths = [x, y]
    halflengths *= 1 // 2
    return CosmoRectangle(center .- halflengths, center .+ halflengths)
end

function CosmoRectangle(center::AbstractVector{<:Number}; halflengths)
    _check_halflengths(center, halflengths)
    return CosmoRectangle(center .- halflengths, center .+ halflengths)
end


"""
    CosmoHypercube(center::AbstractVector{<:Number}, radius::Number)

Return a `CosmoHyperrectangle` centered around `center` whose half side lengths are all equal to `radius`.

The number of dimensions is the length of `center`.

# Examples

```jldoctest
julia> CosmoHypercube([0.0, 0.0, 0.0, 0.0], 1.0)
CosmoHyperrectangle([-1.0, -1.0, -1.0, -1.0], [1.0, 1.0, 1.0, 1.0])
```
"""
function CosmoHypercube(center::AbstractVector{<:Number}, radius::Number)
    return CosmoHyperrectangle(center .- radius, center .+ radius)
end

"""
    CosmoCube(center::AbstractVector{<:Number}, radius::Number)

Return a cubic `CosmoCuboid` centered around `center` whose half side lengths are all equal to `radius`.

# Examples

```jldoctest
julia> CosmoCube([0.0, 0.0, 0.0], 1.0)
CosmoCuboid([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
```
"""
function CosmoCube(center::AbstractVector{<:Number}, radius::Number)
    return CosmoCuboid(center .- radius, center .+ radius)
end

"""
    CosmoSquare(center::AbstractVector{<:Number}, radius::Number)

Return a square `CosmoRectangle` centered around `center` whose half side lengths are both equal to `radius`.

# Examples

```jldoctest
julia> CosmoSquare([0.0, 0.0], 1.0)
CosmoRectangle([-1.0, -1.0], [1.0, 1.0])
```
"""
function CosmoSquare(center::AbstractVector{<:Number}, radius::Number)
    return CosmoRectangle(center .- radius, center .+ radius)
end


for (name, N) in zip([:CosmoHyperrectangle, :CosmoCuboid, :CosmoRectangle], [:(length(lowerleft)), 3, 2])
    quote
        function $name(
            lowerleft::AbstractVector{T1},
            upperright::AbstractVector{T2},
        ) where {T1<:Number,T2<:Number}
            T = promote_type(T1, T2) |> float
            return CosmoHyperrectangle{T,$N}(lowerleft, upperright)
        end
    end |> eval
end

function CosmoHyperrectangle{T}(lowerleft::AbstractVector, upperright::AbstractVector) where {T<:Number}
    return CosmoHyperrectangle{T,length(lowerleft)}(lowerleft, upperright)
end


geometry_enclosing_corners(r::CosmoHyperrectangle) = r.lowerleft, r.upperright

function geometry_enclosing_center(r::CosmoHyperrectangle)
    return 1 // 2 .* (r.lowerleft .+ r.upperright)
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, r::CosmoCuboid{T}) where {T}
    _check_mask_dims(mask, pos, 3)
    Base.require_one_based_indexing(pos)
    @inbounds @views @. mask = (
        (r.lowerleft[1] ≤ pos[1, :] ≤ r.upperright[1]) & (r.lowerleft[2] ≤ pos[2, :] ≤ r.upperright[2]) &
        (r.lowerleft[3] ≤ pos[3, :] ≤ r.upperright[3])
    )
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, r::CosmoRectangle{T}) where {T}
    _check_mask_dims(mask, pos, 2)
    Base.require_one_based_indexing(pos)
    @inbounds @views @. mask =
        ((r.lowerleft[1] ≤ pos[1, :] ≤ r.upperright[1]) & (r.lowerleft[2] ≤ pos[2, :] ≤ r.upperright[2]))
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, r::CosmoHyperrectangle{T,N}) where {T,N}
    _check_mask_dims(mask, pos, N)
    Base.require_one_based_indexing(pos)
    @inbounds @views @. mask = r.lowerleft[1] ≤ pos[1, :] ≤ r.upperright[1]
    for i in 2:N
        mask .&= @inbounds @views @. r.lowerleft[i] ≤ pos[i, :] ≤ r.upperright[i]
    end
    return mask
end

function translate(r::HR, Δx::AbstractVector{<:Number}) where {HR<:CosmoHyperrectangle}
    HR(r.lowerleft .+ Δx, r.upperright .+ Δx)
end


"""
    struct CosmoHypersphere{T<:Number,N} <: AbstractCosmoGeometry
        center::SVector{N,T}
        radius::T
    end

`N`-dimensional hypersphere given by its `center` and `radius`.

If different types are passed to the constructor, the types will be promoted without throwing an error.

# Examples

```jldoctest
julia> s = CosmoHypersphere([0.0, 0.0, 0.0, 0.0], 1.0);

julia> CosmoParticles.mask_in([0.5 2.0; 0.5 0.0; 0.0 0.0; 0.0 0.0], s)
2-element BitVector:
 1
 0
```
"""
struct CosmoHypersphere{T<:Number,N} <: AbstractCosmoGeometry
    center::SVector{N,T}
    radius::T

    function CosmoHypersphere{T,N}(center::AbstractVector{<:Number}, radius::Number) where {T<:Number,N}
        length(center) == N ||
            throw(DimensionMismatch("the center must have length $N, got length $(length(center))"))
        new{T,N}(SVector{N,T}(center), convert(T, radius))
    end
end

"""
    CosmoSphere{T} = CosmoHypersphere{T,3}

Alias for a 3D [`CosmoHypersphere`](@ref).
"""
const CosmoSphere{T} = CosmoHypersphere{T,3} where {T<:Number}

"""
    CosmoCircle{T} = CosmoHypersphere{T,2}

Alias for a 2D [`CosmoHypersphere`](@ref).
"""
const CosmoCircle{T} = CosmoHypersphere{T,2} where {T<:Number}

for (name, N) in zip([:CosmoHypersphere, :CosmoSphere, :CosmoCircle], [:(length(center)), 3, 2])
    quote
        function $name(center::AbstractVector{T1}, radius::T2) where {T1<:Number,T2<:Number}
            T = promote_type(T1, T2) |> float
            return CosmoHypersphere{T,$N}(center, radius)
        end
    end |> eval
end

CosmoHypersphere{T}(center::AbstractVector, radius) where {T<:Number} =
    CosmoHypersphere{T,length(center)}(center, radius)

"""
    CosmoHypersphere(N::Integer, radius::Number)

Returns an `N`-dimensional hypersphere around the origin.

# Examples

```jldoctest
julia> CosmoHypersphere(4, 1.0)
CosmoHypersphere([0.0, 0.0, 0.0, 0.0], 1.0)
```
"""
function CosmoHypersphere(N::Integer, radius::Number)
    return CosmoHypersphere(zeros(typeof(float(radius)), N), radius)
end

"""
    CosmoSphere(radius::Number)

Returns a sphere around the origin.

# Examples

```jldoctest
julia> CosmoSphere(2.0)
CosmoSphere([0.0, 0.0, 0.0], 2.0)
```
"""
function CosmoSphere(radius::Number)
    return CosmoSphere(zeros(typeof(float(radius)), 3), radius)
end

"""
    CosmoCircle(radius::Number)

Returns a circle around the origin.

# Examples

```jldoctest
julia> CosmoCircle(2.0)
CosmoCircle([0.0, 0.0], 2.0)
```
"""
function CosmoCircle(radius::Number)
    return CosmoCircle(zeros(typeof(float(radius)), 2), radius)
end


function Base.:(==)(s1::CosmoHypersphere, s2::CosmoHypersphere)
    return s1.center == s2.center && s1.radius == s2.radius
end

function Base.hash(s::CosmoHypersphere, h::UInt)
    return _hash_normalized(s.radius, _hash_normalized(s.center, hash(:CosmoHypersphere, h)))
end


function geometry_enclosing_corners(s::CosmoHypersphere)
    return s.center .- s.radius, s.center .+ s.radius
end

geometry_enclosing_center(s::CosmoHypersphere) = s.center

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, s::CosmoSphere)
    _check_mask_dims(mask, pos, 3)
    Base.require_one_based_indexing(pos)
    r² = s.radius^2
    if all(==(0), s.center)
        @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 + pos[3, :]^2 ≤ r²
    else
        @inbounds @views @. mask =
            (pos[1, :] - s.center[1])^2 + (pos[2, :] - s.center[2])^2 + (pos[3, :] - s.center[3])^2 ≤ r²
    end
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, circle::CosmoCircle)
    _check_mask_dims(mask, pos, 2)
    Base.require_one_based_indexing(pos)
    r² = circle.radius^2
    if all(==(0), circle.center)
        @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 ≤ r²
    else
        @inbounds @views @. mask = (pos[1, :] - circle.center[1])^2 + (pos[2, :] - circle.center[2])^2 ≤ r²
    end
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, s::CosmoHypersphere{T,N}) where {T,N}
    _check_mask_dims(mask, pos, N)
    Base.require_one_based_indexing(pos)
    r² = s.radius^2
    center = s.center
    if all(==(0), center)
        @inbounds for i in eachindex(mask)
            s = (pos[1, i] - center[1])^2
            for j in 2:length(center)
                s += pos[j, i]^2
            end
            mask[i] = s ≤ r²
        end
    else
        @inbounds for i in eachindex(mask)
            s = (pos[1, i] - center[1])^2
            for j in 2:length(center)
                s += (pos[j, i] - center[j])^2
            end
            mask[i] = s ≤ r²
        end
    end
    return mask
end

function translate(s::HS, Δx::AbstractVector{<:Number}) where {HS<:CosmoHypersphere}
    HS(s.center .+ Δx, s.radius)
end

function rotate(s::HS, rotmat::AbstractMatrix{<:Real}) where {HS<:CosmoHypersphere}
    HS(rotmat * s.center, s.radius)
end


"""
    struct CosmoCylinder{T<:Number} <: AbstractCosmoGeometry
        startpos::SVector{3,T}
        endpos::SVector{3,T}
        radius::T
    end

Cylinder given by its end points `startpos` and `endpos` and `radius`.

If different types are passed to the constructor, the types will be promoted without throwing an error.

# Examples

```jldoctest
julia> c = CosmoCylinder([0.0, 0.0, 0.0], [0.0, 0.0, 2.0], 1.0);

julia> CosmoParticles.mask_in([0.5 0.5; 0.0 0.0; 1.0 3.0], c)
2-element BitVector:
 1
 0
```
"""
struct CosmoCylinder{T<:Number} <: AbstractCosmoGeometry
    startpos::SVector{3,T}
    endpos::SVector{3,T}
    radius::T

    function CosmoCylinder{T}(startpos::AbstractVector, endpos::AbstractVector, radius) where {T<:Number}
        length(startpos) == length(endpos) == 3 || throw(
            DimensionMismatch(
                "the start and end positions must have length 3, got lengths $(length(startpos)) and $(length(endpos))",
            ),
        )
        return new{T}(SVector{3,T}(startpos), SVector{3,T}(endpos), convert(T, radius))
    end
end

function CosmoCylinder(
    startpos::AbstractVector{T1},
    endpos::AbstractVector{T2},
    radius::T3,
) where {T1<:Number,T2<:Number,T3<:Number}
    T = promote_type(T1, T2, T3) |> float
    return CosmoCylinder{T}(startpos, endpos, radius)
end

# Swapping the end points gives the same cylinder.
function Base.:(==)(c1::CosmoCylinder, c2::CosmoCylinder)
    c1.radius == c2.radius || return false
    return c1.startpos == c2.startpos && c1.endpos == c2.endpos ||
           c1.startpos == c2.endpos && c1.endpos == c2.startpos
end


function geometry_enclosing_corners(c::CosmoCylinder)
    # https://iquilezles.org/www/articles/diskbbox/diskbbox.htm
    a = c.endpos .- c.startpos
    @inbounds a² = a[1]^2 + a[2]^2 + a[3]^2
    e = @. c.radius * sqrt(1 - a^2 / a²)
    return (min.(c.startpos .- e, c.endpos .- e), max.(c.startpos .+ e, c.endpos .+ e))
end

geometry_enclosing_center(c::CosmoCylinder) = 1 // 2 .* (c.startpos .+ c.endpos)

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, c::CosmoCylinder)
    # https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
    _check_mask_dims(mask, pos, 3)
    Base.require_one_based_indexing(pos)
    pos1 = @view pos[1, :]
    pos2 = @view pos[2, :]
    pos3 = @view pos[3, :]
    p1 = c.startpos
    p2 = c.endpos
    Δ = p2 .- p1
    norm2Δ = Δ[1]^2 + Δ[2]^2 + Δ[3]^2
    r² = c.radius^2
    @inbounds @. mask =
        ((pos1 - p1[1]) * Δ[1] + (pos2 - p1[2]) * Δ[2] + (pos3 - p1[3]) * Δ[3] ≥ 0) &
        ((pos1 - p2[1]) * Δ[1] + (pos2 - p2[2]) * Δ[2] + (pos3 - p2[3]) * Δ[3] ≤ 0) & (
            (
                ((pos2 - p1[2]) * Δ[3] - (pos3 - p1[3]) * Δ[2])^2 +
                ((pos3 - p1[3]) * Δ[1] - (pos1 - p1[1]) * Δ[3])^2 +
                ((pos1 - p1[1]) * Δ[2] - (pos2 - p1[2]) * Δ[1])^2
            ) / norm2Δ ≤ r²
        )
end

function translate(c::C, Δx::AbstractVector{<:Number}) where {C<:CosmoCylinder}
    C(c.startpos .+ Δx, c.endpos .+ Δx, c.radius)
end

function rotate(c::C, rotmat::AbstractMatrix{<:Real}) where {C<:CosmoCylinder}
    C(rotmat * c.startpos, rotmat * c.endpos, c.radius)
end


"""
    struct CosmoStandingCylinder{T<:Number} <: AbstractCosmoGeometry
        center::SVector{3,T}
        height::T
        radius::T
    end

Standing cylinder given by its `center`, `height`, and `radius`.

Standing means that the cylinder is oriented such that its axis is aligned with the z axis.
If different types are passed to the constructor, the types will be promoted without throwing an error.

# Examples

```jldoctest
julia> c = CosmoStandingCylinder([0.0, 0.0, 0.0], 2.0, 1.0);

julia> CosmoParticles.mask_in([0.5 0.5; 0.0 0.0; 0.5 3.0], c)
2-element BitVector:
 1
 0
```
"""
struct CosmoStandingCylinder{T<:Number} <: AbstractCosmoGeometry
    center::SVector{3,T}
    height::T
    radius::T

    function CosmoStandingCylinder{T}(center::AbstractVector, height, radius) where {T<:Number}
        length(center) == 3 ||
            throw(DimensionMismatch("the center must have length 3, got length $(length(center))"))
        return new{T}(SVector{3,T}(center), convert(T, height), convert(T, radius))
    end
end

function CosmoStandingCylinder(
    center::AbstractVector{T1},
    height::T2,
    radius::T3,
) where {T1<:Number,T2<:Number,T3<:Number}
    T = promote_type(T1, T2, T3) |> float
    return CosmoStandingCylinder{T}(center, height, radius)
end

function Base.:(==)(c1::CosmoStandingCylinder, c2::CosmoStandingCylinder)
    return c1.center == c2.center && c1.height == c2.height && c1.radius == c2.radius
end


function geometry_enclosing_corners(c::CosmoStandingCylinder)
    halfsize = SVector(c.radius, c.radius, 1 // 2 * c.height)
    return c.center .- halfsize, c.center .+ halfsize
end

geometry_enclosing_center(c::CosmoStandingCylinder) = c.center

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, c::CosmoStandingCylinder)
    _check_mask_dims(mask, pos, 3)
    Base.require_one_based_indexing(pos)
    r² = c.radius^2
    halfheight = 1 // 2 * c.height
    @inbounds @views @. mask =
        ((pos[1, :] - c.center[1])^2 + (pos[2, :] - c.center[2])^2 ≤ r²) &
        (c.center[3] - halfheight ≤ pos[3, :] ≤ c.center[3] + halfheight)
end

function translate(c::C, Δx::AbstractVector{<:Number}) where {C<:CosmoStandingCylinder}
    C(c.center .+ Δx, c.height, c.radius)
end

function rotate(c::C, rotmat::AbstractMatrix{<:Real}) where {C<:CosmoStandingCylinder}
    rotate(CosmoCylinder(c), rotmat)
end


"""
    CosmoCylinder(c::CosmoStandingCylinder)

Create a cylinder from a standing cylinder.

# Examples

```jldoctest
julia> CosmoCylinder(CosmoStandingCylinder([0.0, 0.0, 0.0], 2.0, 1.0))
CosmoCylinder([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 1.0)
```
"""
function CosmoCylinder(c::CosmoStandingCylinder)
    halfheight = 1 // 2 * c.height
    startpos = [c.center[1], c.center[2], c.center[3] - halfheight]
    endpos = [c.center[1], c.center[2], c.center[3] + halfheight]
    return CosmoCylinder(startpos, endpos, c.radius)
end

"""
    CosmoStandingCylinder(c::CosmoCylinder)

Create a standing cylinder from a cylinder.

The end positions of the cylinder have to have the same x and y coordinates.

# Examples

```jldoctest
julia> c = CosmoCylinder([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 1.0);

julia> CosmoStandingCylinder(c)
CosmoStandingCylinder([0.0, 0.0, 0.0], 2.0, 1.0)
```
"""
function CosmoStandingCylinder(c::CosmoCylinder)
    c.startpos[1] == c.endpos[1] && c.startpos[2] == c.endpos[2] || throw(
        ArgumentError(
            "the cylinder must be aligned with the z axis, got start position $(c.startpos) and end position $(c.endpos)",
        ),
    )
    center = geometry_enclosing_center(c)
    height = abs(c.endpos[3] - c.startpos[3])
    return CosmoStandingCylinder(center, height, c.radius)
end

function Base.:(==)(c1::CosmoCylinder, c2::CosmoStandingCylinder)
    return geometry_enclosing_center(c1) == geometry_enclosing_center(c2) &&
           abs.(c1.startpos[3] - c1.endpos[3]) == c2.height &&
           c1.radius == c2.radius &&
           @views c1.startpos[1:2] == c1.endpos[1:2] # c1 is actually standing
end
Base.:(==)(c1::CosmoStandingCylinder, c2::CosmoCylinder) = ==(c2, c1)

# A cylinder along the z axis can equal a standing cylinder, so both hash its center, height, and radius.
function Base.hash(c::CosmoCylinder, h::UInt)
    if @views c.startpos[1:2] == c.endpos[1:2]
        return _hash_standing_cylinder(
            geometry_enclosing_center(c),
            abs(c.endpos[3] - c.startpos[3]),
            c.radius,
            h,
        )
    end
    # The end points are hashed in a fixed order, since swapping them gives the same cylinder.
    p1, p2 = c.startpos .+ zero.(c.startpos), c.endpos .+ zero.(c.endpos)
    if isless(p2, p1)
        p1, p2 = p2, p1
    end
    return _hash_normalized(c.radius, hash(p2, hash(p1, hash(:CosmoCylinder, h))))
end
Base.hash(c::CosmoStandingCylinder, h::UInt) = _hash_standing_cylinder(c.center, c.height, c.radius, h)

function _hash_standing_cylinder(center, height, radius, h::UInt)
    h = _hash_normalized(center, hash(:CosmoStandingCylinder, h))
    return _hash_normalized(radius, _hash_normalized(height, h))
end

function Base.isapprox(c1::CosmoCylinder, c2::CosmoStandingCylinder; kwargs...)
    return isapprox(geometry_enclosing_center(c1), geometry_enclosing_center(c2); kwargs...) &&
           isapprox(abs.(c1.startpos[3] - c1.endpos[3]), c2.height; kwargs...) &&
           isapprox(c1.radius, c2.radius; kwargs...) &&
           @views isapprox(c1.startpos[1:2], c1.endpos[1:2]; kwargs...) # c1 is actually standing
end
Base.isapprox(c1::CosmoStandingCylinder, c2::CosmoCylinder; kwargs...) = isapprox(c2, c1; kwargs...)

for G in (CosmoUnionGeometry, CosmoIntersectGeometry, CosmoDiffGeometry, Rotated, Translated)
    @eval Base.:(==)(g1::$G, g2::$G) = _fields_equal(g1, g2)
    @eval Base.hash(g::$G, h::UInt) = _fields_hash(g, h)
end

# Geometries are shown as the call that constructs them.
function _show_call(io::IO, f, args...; kwargs...)
    print(io, f, "(")
    join(io, (_repr_arg(io, arg) for arg in args), ", ")
    if !isempty(kwargs)
        print(io, "; ")
        join(io, ("$k=$(_repr_arg(io, v))" for (k, v) in kwargs), ", ")
    end
    print(io, ")")
    return nothing
end

# Quantities are shown with the string macro of Unitful, e.g. as `[1.0, 2.0]u"m"`.
_repr_arg(io::IO, x) = repr(x; context=io)
_repr_arg(io::IO, x::Quantity) = string(repr(ustrip(x); context=io), _unit_string(unit(x)))
function _repr_arg(io::IO, x::AbstractArray{<:Quantity})
    isconcretetype(eltype(x)) || return repr(x; context=io)
    return string(repr(ustrip.(x); context=io), _unit_string(unit(eltype(x))))
end
# Unitful writes exponents as superscripts on macOS by default, which the string macro cannot parse.
_unit_string(u) = string("u\"", replace(sprint(show, u; context=:fancy_exponent => false), " " => "*"), "\"")

function Base.show(io::IO, r::CosmoHyperrectangle{T,N}) where {T,N}
    name = N == 3 ? "CosmoCuboid" : (N == 2 ? "CosmoRectangle" : "CosmoHyperrectangle")
    _show_call(io, name, r.lowerleft, r.upperright)
end
function Base.show(io::IO, s::CosmoHypersphere{T,N}) where {T,N}
    name = N == 3 ? "CosmoSphere" : (N == 2 ? "CosmoCircle" : "CosmoHypersphere")
    _show_call(io, name, s.center, s.radius)
end
Base.show(io::IO, c::CosmoCylinder) = _show_call(io, "CosmoCylinder", c.startpos, c.endpos, c.radius)
Base.show(io::IO, c::CosmoStandingCylinder) =
    _show_call(io, "CosmoStandingCylinder", c.center, c.height, c.radius)
Base.show(io::IO, g::CosmoUnionGeometry) = _show_call(io, "union", g.geos...)
Base.show(io::IO, g::CosmoIntersectGeometry) = _show_call(io, "intersect", g.geos...)
Base.show(io::IO, g::CosmoDiffGeometry) = _show_call(io, "setdiff", g.geo, g.geos...)
Base.show(io::IO, g::Rotated) = _show_call(io, "rotate", g.geo, g.rotmat)
Base.show(io::IO, g::Translated) = _show_call(io, "translate", g.geo, g.Δx)
