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
"""
abstract type AbstractCosmoGeometry end

"""
    geometry_enclosing_corners(geo::AbstractCosmoGeometry)

Return the lower left and upper right corners of the enclosing box of the geometry as a tuple of vectors.

The enclosing box is not necessarily the tightest fitting box.

This is not exported.
"""
function geometry_enclosing_corners(::AbstractCosmoGeometry) end

"""
    geometry_enclosing_center(geo::AbstractCosmoGeometry)

Return the center of the enclosing box of the geometry as a vector.

This is not exported.
"""
function geometry_enclosing_center(geo::AbstractCosmoGeometry)
    corners = geometry_enclosing_corners(geo)
    return 1 // 2 .* (corners[1] .+ corners[2])
end

@doc raw"""
    mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)

Return the `BitVector` mask of the positions (``\mathrm{dims} × N``) located within the geometry.

This is not exported.
"""
function mask_in!(::BitVector, ::AbstractMatrix{<:Number}, ::AbstractCosmoGeometry) end

@doc raw"""
    mask_in(pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)

Return the `BitArray` mask of the positions (``\mathrm{dims} × N``) located within the geometry.

Calls [`mask_in!`](@ref) by default.

This is not exported.
"""
function mask_in(pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)
    mask = BitVector(undef, size(pos, 2))
    return mask_in!(mask, pos, geo)
end



"""
    struct CosmoUnionGeometry <: AbstractCosmoGeometry
        geos::Tuple
    end

Union of multiple geometries.

This union geometry represents all of the volumes contained within the geometries. 

This is not exported.
"""
struct CosmoUnionGeometry <: AbstractCosmoGeometry
    geos::Tuple

    function CosmoUnionGeometry(geos)
        new(Tuple(geos))
    end
end

"""
    CosmoUnionGeometry(geos::AbstractVector) 
    CosmoUnionGeometry(geos::AbstractCosmoGeometry...) 

Create a union of geometries from the passed geometries.

This is not exported.
"""
CosmoUnionGeometry(geos::AbstractCosmoGeometry...) = CosmoUnionGeometry(geos)

"""
    Base.union(geos::AbstractCosmoGeometry...) 

Create a union of geometries from the passed geometries.

# Examples
```julia
using CosmoParticles
using CosmoParticles: mask_in

pos::Matrix # 3×n

g1 = CosmoSphere([1, 2, 3], 4)
g2 = CosmoSphere([3, 2, 3], 4)
geo = union(g1, g2)

mask_in(pos, geo) == mask_in(pos, g1) .| mask_in(pos, g2)
```
"""
Base.union(geos::AbstractCosmoGeometry...) = CosmoUnionGeometry(geos)

function geometry_enclosing_corners(g::CosmoUnionGeometry)
    corners = [geometry_enclosing_corners(geo) for geo in g.geos]
    ndims = length(corners[1][1])
    lowerleft = [minimum(ll[i] for (ll, ur) in corners) for i in 1:ndims]
    upperright = [maximum(ur[i] for (ll, ur) in corners) for i in 1:ndims]

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

This is not exported.
"""
struct CosmoIntersectGeometry <: AbstractCosmoGeometry
    geos::Tuple

    function CosmoIntersectGeometry(geos)
        new(Tuple(geos))
    end
end

"""
    CosmoIntersectGeometry(geos::AbstractVector) 
    CosmoIntersectGeometry(geos::AbstractCosmoGeometry...) 

Create an intersect of geometries from the passed geometries.

This is not exported.
"""
CosmoIntersectGeometry(geos::AbstractCosmoGeometry...) = CosmoIntersectGeometry(geos)

"""
    Base.intersect(geos::AbstractCosmoGeometry...) 

Create an intersect of geometries from the passed geometries.

# Examples
```julia
using CosmoParticles
using CosmoParticles: mask_in

pos::Matrix # 3×n

g1 = CosmoSphere([1, 2, 3], 4)
g2 = CosmoSphere([3, 2, 3], 4)
geo = intersect(g1, g2)

mask_in(pos, geo) == mask_in(pos, g1) .& mask_in(pos, g2)
```
"""
Base.intersect(geos::AbstractCosmoGeometry...) = CosmoIntersectGeometry(geos)

function geometry_enclosing_corners(g::CosmoIntersectGeometry)
    corners = [geometry_enclosing_corners(geo) for geo in g.geos]
    ndims = length(corners[1][1])
    lowerleft = [maximum(ll[i] for (ll, ur) in corners) for i in 1:ndims]
    upperright = [minimum(ur[i] for (ll, ur) in corners) for i in 1:ndims]

    if any(lowerleft .≥ upperright)
        return zeros(ndims), zeros(ndims)
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

This is not exported.
"""
struct CosmoDiffGeometry <: AbstractCosmoGeometry
    geo::AbstractCosmoGeometry
    geos::Tuple

    function CosmoDiffGeometry(geo::AbstractCosmoGeometry, geos)
        new(geo, Tuple(geos))
    end
end

"""
    CosmoDiffGeometry(geo::AbstractCosmoGeometry, geos::AbstractVector) 
    CosmoDiffGeometry(geo::AbstractCosmoGeometry, geos::AbstractCosmoGeometry...) 

Create a difference of geometries from the passed geometries.

This is not exported.
"""
CosmoDiffGeometry(geo::AbstractCosmoGeometry, geos::AbstractCosmoGeometry...) = CosmoDiffGeometry(geo, geos)

"""
    Base.setdiff(geo::AbstractCosmoGeometry, geos::AbstractCosmoGeometry...) 

Create a difference of geometries from the passed geometries.

# Examples
```julia
using CosmoParticles
using CosmoParticles: mask_in

pos::Matrix # 3×n

g1 = CosmoSphere([1, 2, 3], 4)
g2 = CosmoSphere([3, 2, 3], 4)
geo = setdiff(g1, g2)

mask_in(pos, geo) == mask_in(pos, g1) .& .~mask_in(pos, g2)
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
    struct Rotated{G,R} <: AbstractCosmoGeometry where {G<:AbstractCosmoGeometry,R<:AbstractMatrix}
        geo::G
        rotmat::R
    end

Rotated geometry defined by a rotation matrix.

This can be created by [`rotate`](@ref) applied to a geometry.

This is not exported.
"""
struct Rotated{G,R} <: AbstractCosmoGeometry where {G<:AbstractCosmoGeometry,R<:AbstractMatrix}
    geo::G
    rotmat::R
end

"""
    rotation_matrix(rot::Rotated)

Returns the rotation matrix that rotates the original geometry into the rotated state.

This is not exported.
"""
rotation_matrix(rot::Rotated) = rot.rotmat

"""
    rotation_matrix_inv(rot::Rotated)

Returns the rotation matrix that rotates the rotated geometry back to its original orientation.

Use this to rotate objects in the rotated frame of reference to the original frame of reference
of the geometry.

This is not exported.
"""
rotation_matrix_inv(rot::Rotated) = transpose(rot.rotmat)

function geometry_enclosing_corners(rot::Rotated)
    bounding_corners = geometry_enclosing_corners(rot.geo)
    n = length(bounding_corners[1])

    # create matrix of positions of the hypercuboid corners
    pos = Matrix{eltype(bounding_corners[1])}(undef, n, 2^n)
    @inbounds for bitval in axes(pos, 2)
        for i in 0:n-1
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
    posrot = pos * rot.rotmat

    # extract the extrema along each dimension as the lower left and upper right corners
    return minimum(posrot; dims=2), maximum(posrot; dims=2)
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, rot::Rotated)
    mask_in!(mask, matrix_rotate(pos, rotation_matrix_inv(rot)), rot.geo)
end

"""
    rotate(geo::AbstractCosmoGeometry, rotmat::AbstractMatrix{<:Real})

Returns a rotated geometry rotated by the rotation matrix.

Individual geometries may specify that this returns a different type than [`Rotated`](@ref).
"""
rotate(geo::AbstractCosmoGeometry, rotmat::AbstractMatrix{<:Real}) = Rotated(geo, rotmat)
rotate(rot::Rotated, rotmat::AbstractMatrix{<:Real}) = Rotated(rot.geo, rotmat * rot.rotmat)



"""
    struct Translated{G,L} <: AbstractCosmoGeometry where {G<:AbstractCosmoGeometry,L<:AbstractVector}
        geo::G
        Δx::L
    end

Translated geometry defined by a translation vector.

This can be created by [`translate`](@ref) applied to a geometry.

This is not exported.
"""
struct Translated{G,L} <: AbstractCosmoGeometry where {G<:AbstractCosmoGeometry,L<:AbstractVector}
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
"""
translate(geo::AbstractCosmoGeometry, Δx::AbstractVector{<:Number}) = Translated(geo, Δx)
translate(tra::Translated, Δx::AbstractVector{<:Number}) = Translated(tra.geo, tra.Δx .+ Δx)





"""
    struct CosmoHyperrectangle{T,N} <: AbstractCosmoGeometry where {T<:Number}
        lowerleft::Vector{T}
        upperright::Vector{T}
    end

`N`-dimensional hyperrectangle aligned with the coordinate system axes given by its lower left and upper right corners.

The dimensions of space are given by `N`.
"""
struct CosmoHyperrectangle{T,N} <: AbstractCosmoGeometry where {T<:Number}
    lowerleft::Vector{T}
    upperright::Vector{T}

    function CosmoHyperrectangle{T,N}(
        lowerleft::AbstractVector{<:Number},
        upperright::AbstractVector{<:Number},
    ) where {T<:Number,N}
        @assert length(lowerleft) == length(upperright) == N
        new{T,N}(convert(Vector{T}, lowerleft), convert(Vector{T}, upperright))
    end
end

"""
    CosmoHyperrectangle(center::AbstractVector{<:Number}, sidelengths::Number...)

Return a `CosmoHyperrectangle` centered around `center`, with the given `sidelengths`.

The side lengths are given in order of the axes, i.e., x, y or x, y, z for 2D and 3D, respectively.
"""
function CosmoHyperrectangle(center::AbstractVector{<:Number}, sidelengths::Number...)
    halflengths = [sidelengths...]
    halflengths *= 1 // 2
    return CosmoHyperrectangle(center .- halflengths, center .+ halflengths)
end

function Base.:(==)(r1::CosmoHyperrectangle, r2::CosmoHyperrectangle)
    return r1.lowerleft == r2.lowerleft && r1.upperright == r2.upperright
end


"""
    CosmoCuboid{T} = CosmoHyperrectangle{T,3}

Alias for a 3D [`CosmoHyperrectangle`](@ref).
"""
const CosmoCuboid{T} = CosmoHyperrectangle{T,3} where {T<:Number}

function CosmoCuboid(center::AbstractVector{<:Number}, x::Number, y::Number, z::Number)
    halflengths = [x, y, z]
    halflengths *= 1 // 2
    return CosmoCuboid(center .- halflengths, center .+ halflengths)
end


"""
    CosmoRectangle{T} = CosmoHyperrectangle{T,2}

Alias for a 2D [`CosmoHyperrectangle`](@ref).
"""
const CosmoRectangle{T} = CosmoHyperrectangle{T,2} where {T<:Number}

function CosmoRectangle(center::AbstractVector{<:Number}, x::Number, y::Number)
    halflengths = [x, y]
    halflengths *= 1 // 2
    return CosmoRectangle(center .- halflengths, center .+ halflengths)
end


"""
    CosmoHypercube(center::AbstractVector{<:Number}, radius::Number)

Return a cubic `CosmoCuboid` centered around `center`, with equal sidelengths.
"""
function CosmoHypercube(center::AbstractVector{<:Number}, radius::Number)
    return CosmoHyperrectangle(center .- radius, center .+ radius)
end

"""
    CosmoCube(center::AbstractVector{<:Number}, radius::Number)

Return a cubic `CosmoCuboid` centered around `center`, with equal sidelengths.
"""
function CosmoCube(center::AbstractVector{<:Number}, radius::Number)
    return CosmoCuboid(center .- radius, center .+ radius)
end

"""
    CosmoSquare(center::AbstractVector{<:Number}, radius::Number)

Return a square `CosmoRectangle` centered around `center`, with equal sidelengths.
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


geometry_enclosing_corners(r::CosmoHyperrectangle) = r.lowerleft, r.upperright

function geometry_enclosing_center(r::CosmoHyperrectangle)
    return 1 // 2 .* (r.lowerleft .+ r.upperright)
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, r::CosmoCuboid{T}) where {T}
    @assert size(pos, 1) == 3
    @assert size(pos, 2) == length(mask)
    @inbounds @views @. mask = (
        (r.lowerleft[1] ≤ pos[1, :] ≤ r.upperright[1]) &
        (r.lowerleft[2] ≤ pos[2, :] ≤ r.upperright[2]) &
        (r.lowerleft[3] ≤ pos[3, :] ≤ r.upperright[3])
    )
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, r::CosmoRectangle{T}) where {T}
    @assert size(pos, 1) == 2
    @assert size(pos, 2) == length(mask)
    @inbounds @views @. mask =
        ((r.lowerleft[1] ≤ pos[1, :] ≤ r.upperright[1]) & (r.lowerleft[2] ≤ pos[2, :] ≤ r.upperright[2]))
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, r::CosmoHyperrectangle{T,N}) where {T,N}
    @assert size(pos, 1) == N
    @assert size(pos, 2) == length(mask)
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
    struct CosmoHypersphere{T,N} <: AbstractCosmoGeometry where {T<:Number}
        center::Vector{T}
        radius::T
    end

`N`-dimensional hypersphere given by its `center` and `radius`.

If different types are passed to the constructor, the types will be promoted without throwing an error.
"""
struct CosmoHypersphere{T,N} <: AbstractCosmoGeometry where {T<:Number}
    center::Vector{T}
    radius::T

    function CosmoHypersphere{T,N}(center::AbstractVector{<:Number}, radius::Number) where {T<:Number,N}
        @assert length(center) == N
        new{T,N}(convert(Vector{T}, center), convert(T, radius))
    end
end

"""
    CosmoSphere{T} = CosmoHypersphere{T,3}

Alias for a 2D [`CosmoHypersphere`](@ref).
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

"""
    CosmoHypersphere(N::Integer, radius::T)

Returns an `N`-dimensional hypersphere around the origin.
"""
function CosmoHypersphere(N::Integer, radius::T) where {T<:Number}
    Tfloat = float(T)
    return CosmoHypersphere{Tfloat,N}(zeros(Tfloat, N), radius)
end

"""
    CosmoSphere(radius::T)

Returns a sphere around the origin.
"""
function CosmoSphere(radius::T) where {T<:Number}
    Tfloat = float(T)
    return CosmoHypersphere{Tfloat,3}(zeros(Tfloat, 3), radius)
end

"""
    CosmoCircle(radius::T)

Returns a circle around the origin.
"""
function CosmoCircle(radius::T) where {T<:Number}
    Tfloat = float(T)
    return CosmoHypersphere{Tfloat,2}(zeros(Tfloat, 2), radius)
end


function Base.:(==)(s1::CosmoHypersphere, s2::CosmoHypersphere)
    return s1.center == s2.center && s1.radius == s2.radius
end


function geometry_enclosing_corners(s::CosmoHypersphere)
    return s.center .- s.radius, s.center .+ s.radius
end

geometry_enclosing_center(s::CosmoHypersphere) = s.center

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, s::CosmoSphere)
    @assert size(pos, 1) == 3
    @assert size(pos, 2) == length(mask)
    r² = s.radius^2
    if all(==(0), s.center)
        @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 + pos[3, :]^2 ≤ r²
    else
        @inbounds @views @. mask =
            (pos[1, :] - s.center[1])^2 + (pos[2, :] - s.center[2])^2 + (pos[3, :] - s.center[3])^2 ≤ r²
    end
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, circle::CosmoCircle)
    @assert size(pos, 1) == 2
    @assert size(pos, 2) == length(mask)
    r² = circle.radius^2
    if all(==(0), circle.center)
        @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 ≤ r²
    else
        @inbounds @views @. mask = (pos[1, :] - circle.center[1])^2 + (pos[2, :] - circle.center[2])^2 ≤ r²
    end
end

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, s::CosmoHypersphere{T,N}) where {T,N}
    @assert size(pos, 1) == N
    @assert size(pos, 2) == length(mask)
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
    struct CosmoCylinder{T} <: AbstractCosmoGeometry where {T<:Number}
        startpos::Vector{T}
        endpos::Vector{T}
        radius::T
    end

Cylinder given by its end points `startpos` and `endpos` and `radius`.

If different types are passed to the constructor, the types will be promoted without throwing an error.
"""
struct CosmoCylinder{T} <: AbstractCosmoGeometry where {T<:Number}
    startpos::Vector{T}
    endpos::Vector{T}
    radius::T

    function CosmoCylinder{T}(startpos::AbstractVector, endpos::AbstractVector, radius) where {T<:Number}
        @assert length(startpos) == 3 && length(endpos) == 3
        return new{T}(convert(Vector{T}, startpos), convert(Vector{T}, endpos), convert(T, radius))
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

function Base.:(==)(c1::CosmoCylinder, c2::CosmoCylinder)
    return c1.startpos == c2.startpos && c1.endpos == c2.endpos && c1.radius == c2.radius
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
    @assert size(pos, 1) == 3
    @assert size(pos, 2) == length(mask)
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
        ((pos1 - p2[1]) * Δ[1] + (pos2 - p2[2]) * Δ[2] + (pos3 - p2[3]) * Δ[3] ≤ 0) &
        (
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
    struct CosmoStandingCylinder{T} <: AbstractCosmoGeometry where {T<:Number}
        center::Vector{T}
        height::T
        radius::T
    end

Standing cylinder given by its `center`, `height`, and `radius`.

Standing means that the cylinder is oriented such that its axis is aligned with the z axis.
If different types are passed to the constructor, the types will be promoted without throwing an error.
"""
struct CosmoStandingCylinder{T} <: AbstractCosmoGeometry where {T<:Number}
    center::Vector{T}
    height::T
    radius::T

    function CosmoStandingCylinder{T}(center::AbstractVector, height, radius) where {T<:Number}
        @assert length(center) == 3
        return new{T}(convert(Vector{T}, center), convert(T, height), convert(T, radius))
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
    c.center .- [c.radius, c.radius, 1 // 2 * c.height], c.center .+ [c.radius, c.radius, 1 // 2 * c.height]
end

geometry_enclosing_center(c::CosmoStandingCylinder) = c.center

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, c::CosmoStandingCylinder)
    @assert size(pos, 1) == 3
    @assert size(pos, 2) == length(mask)
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
"""
function CosmoStandingCylinder(c::CosmoCylinder)
    @assert c.startpos[1] == c.endpos[1] && c.startpos[2] == c.endpos[2]
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

function Base.isapprox(c1::CosmoCylinder, c2::CosmoStandingCylinder; kwargs...)
    return isapprox(geometry_enclosing_center(c1), geometry_enclosing_center(c2); kwargs...) &&
           isapprox(abs.(c1.startpos[3] - c1.endpos[3]), c2.height; kwargs...) &&
           isapprox(c1.radius, c2.radius; kwargs...) &&
           @views isapprox(c1.startpos[1:2], c1.endpos[1:2]; kwargs...) # c1 is actually standing
end
Base.isapprox(c1::CosmoStandingCylinder, c2::CosmoCylinder; kwargs...) = isapprox(c2, c1; kwargs...)
