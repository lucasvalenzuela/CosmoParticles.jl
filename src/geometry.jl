"""
    abstract type AbstractCosmoGeometry end

Abstract type for (multi-dimensional) geometry volumes, particularly for filtering with [`filter`](@ref).

Any subtypes of `AbstractCosmoGeometry` have to implement the following methods:
- [`CosmoParticles.geometry_enclosing_corners`](@ref)
- [`CosmoParticles.geometry_enclosing_center`](@ref)
- [`CosmoParticles.mask_in`](@ref)
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
function geometry_enclosing_center(::AbstractCosmoGeometry) end

@doc raw"""
    mask_in(pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)

Return the `BitArray` mask of the positions (``\mathrm{dims} × N``) located within the geometry.

This is not exported.
"""
function mask_in(::AbstractMatrix{<:Number}, ::AbstractCosmoGeometry) end


"""
    struct CosmoHyperrectangle{T,N} <: AbstractCosmoGeometry where {T<:Number}
        lowerleft::Vector{T}
        upperright::Vector{T}
    end

Hyperrectangle aligned with the coordinate system axes given by its lower left and upper right corners.

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
            T = promote_type(T1, T2)
            return CosmoHyperrectangle{T,$N}(lowerleft, upperright)
        end
    end |> eval
end


geometry_enclosing_corners(r::CosmoHyperrectangle) = r.lowerleft, r.upperright

function geometry_enclosing_center(r::CosmoHyperrectangle)
    return 1 // 2 .* (r.lowerleft .+ r.upperright)
end

function mask_in(pos::AbstractMatrix{<:Number}, r::CosmoCuboid{T}) where {T}
    @assert size(pos, 1) == 3
    @inbounds @views @. (
        (r.lowerleft[1] ≤ pos[1, :] ≤ r.upperright[1]) &
        (r.lowerleft[2] ≤ pos[2, :] ≤ r.upperright[2]) &
        (r.lowerleft[3] ≤ pos[3, :] ≤ r.upperright[3])
    )
end

function mask_in(pos::AbstractMatrix{<:Number}, r::CosmoRectangle{T}) where {T}
    @assert size(pos, 1) == 2
    @inbounds @views @. (
        (r.lowerleft[1] ≤ pos[1, :] ≤ r.upperright[1]) & (r.lowerleft[2] ≤ pos[2, :] ≤ r.upperright[2])
    )
end

function mask_in(pos::AbstractMatrix{<:Number}, r::CosmoHyperrectangle{T,N}) where {T,N}
    @assert size(pos, 1) == N
    mask = @inbounds @views @. r.lowerleft[1] ≤ pos[1, :] ≤ r.upperright[1]
    for i in 2:N
        mask .&= @inbounds @views @. r.lowerleft[i] ≤ pos[i, :] ≤ r.upperright[i]
    end
    return mask
end


"""
    struct CosmoHypersphere{T,N} <: AbstractCosmoGeometry where {T<:Number}
        center::Vector{T}
        radius::T
    end

Hypersphere given by its `center` and `radius`.

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
            T = promote_type(T1, T2)
            return CosmoHypersphere{T,$N}(center, radius)
        end
    end |> eval
end


function Base.:(==)(s1::CosmoHypersphere, s2::CosmoHypersphere)
    return s1.center == s2.center && s1.radius == s2.radius
end


function geometry_enclosing_corners(s::CosmoHypersphere)
    return s.center .- s.radius, s.center .+ s.radius
end

geometry_enclosing_center(s::CosmoHypersphere) = s.center

function mask_in(pos::AbstractMatrix{<:Number}, s::CosmoSphere)
    @assert size(pos, 1) == 3
    r² = s.radius^2
    if all(==(0), s.center)
        @inbounds @views @. pos[1, :]^2 + pos[2, :]^2 + pos[3, :]^2 ≤ r²
    else
        @inbounds @views @. (pos[1, :] - s.center[1])^2 +
                            (pos[2, :] - s.center[2])^2 +
                            (pos[3, :] - s.center[3])^2 ≤ r²
    end
end

function mask_in(pos::AbstractMatrix{<:Number}, circle::CosmoCircle)
    @assert size(pos, 1) == 2
    r² = circle.radius^2
    if all(==(0), circle.center)
        @inbounds @views @. pos[1, :]^2 + pos[2, :]^2 ≤ r²
    else
        @inbounds @views @. (pos[1, :] - circle.center[1])^2 + (pos[2, :] - circle.center[2])^2 ≤ r²
    end
end

function mask_in(pos::AbstractMatrix{<:Number}, s::CosmoHypersphere{T,N}) where {T,N}
    @assert size(pos, 1) == N
    r² = s.radius^2
    center = s.center
    mask = BitArray(undef, size(pos, 2))
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
    T = promote_type(T1, T2, T3)
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

function mask_in(pos::AbstractMatrix{<:Number}, c::CosmoCylinder)
    # https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
    @assert size(pos, 1) == 3
    pos1 = @view pos[1, :]
    pos2 = @view pos[2, :]
    pos3 = @view pos[3, :]
    p1 = c.startpos
    p2 = c.endpos
    Δ = p2 .- p1
    norm2Δ = Δ[1]^2 + Δ[2]^2 + Δ[3]^2
    r² = c.radius^2
    @inbounds @. ((pos1 - p1[1]) * Δ[1] + (pos2 - p1[2]) * Δ[2] + (pos3 - p1[3]) * Δ[3] ≥ 0) &
                 ((pos1 - p2[1]) * Δ[1] + (pos2 - p2[2]) * Δ[2] + (pos3 - p2[3]) * Δ[3] ≤ 0) &
                 (
                     (
                         ((pos2 - p1[2]) * Δ[3] - (pos3 - p1[3]) * Δ[2])^2 +
                         ((pos3 - p1[3]) * Δ[1] - (pos1 - p1[1]) * Δ[3])^2 +
                         ((pos1 - p1[1]) * Δ[2] - (pos2 - p1[2]) * Δ[1])^2
                     ) / norm2Δ ≤ r²
                 )
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
    T = promote_type(T1, T2, T3)
    return CosmoStandingCylinder{T}(center, height, radius)
end

function Base.:(==)(c1::CosmoStandingCylinder, c2::CosmoStandingCylinder)
    return c1.center == c2.center && c1.height == c2.height && c1.radius == c2.radius
end


function geometry_enclosing_corners(c::CosmoStandingCylinder)
    c.center .- [c.radius, c.radius, 1 // 2 * c.height], c.center .+ [c.radius, c.radius, 1 // 2 * c.height]
end

geometry_enclosing_center(c::CosmoStandingCylinder) = c.center

function mask_in(pos::AbstractMatrix{<:Number}, c::CosmoStandingCylinder)
    @assert size(pos, 1) == 3
    r² = c.radius^2
    halfheight = 1 // 2 * c.height
    @inbounds @views @. ((pos[1, :] - c.center[1])^2 + (pos[2, :] - c.center[2])^2 ≤ r²) &
                        (c.center[3] - halfheight ≤ pos[3, :] ≤ c.center[3] + halfheight)
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
