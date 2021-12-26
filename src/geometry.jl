abstract type AbstractCosmoGeometry end

"""
    geometry_enclosing_corners(geo::AbstractCosmoGeometry)

Return the lower left and upper right corners of the enclosing box of the geometry.
"""
function geometry_enclosing_corners(::AbstractCosmoGeometry) end

"""
    geometry_enclosing_center(geo::AbstractCosmoGeometry)

Return the center of the enclosing box of the geometry.
"""
function geometry_enclosing_center(::AbstractCosmoGeometry) end

@doc raw"""
    mask_in(pos::AbstractMatrix{<:Number}, geo::AbstractCosmoGeometry)

Return the `BitArray` mask of the positions (``\mathrm{dims} \times N``) located within the geometry.
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

    function CosmoHyperrectangle{T,N}(lowerleft::Vector{T}, upperright::Vector{T}) where {T<:Number,N}
        @assert length(lowerleft) == length(upperright) == N
        new{T,N}(lowerleft, upperright)
    end
end

"""
    CosmoHyperrectangle(center::AbstractVector{<:Number}, sidelengths::Number...)

Return a `CosmoHyperrectangle` centered around `center`, with the given `sidelengths`.

The side lengths are given in order of the axes, i.e., x, y or x, y, z for 2D and 3D, respectively.
"""
function CosmoHyperrectangle(center::Vector{<:Number}, sidelengths::Number...)
    halflengths = [sidelengths...]
    halflengths *= 1 // 2
    return CosmoHyperrectangle(center .- halflengths, center .+ halflengths)
end

function Base.:(==)(r1::CosmoHyperrectangle, r2::CosmoHyperrectangle)
    return r1.lowerleft == r2.lowerleft && r1.upperright == r2.upperright
end


"""
    CosmoCuboid

Alias for a 3D [`CosmoHyperrectangle`](@ref).
"""
const CosmoCuboid{T} = CosmoHyperrectangle{T,3} where {T<:Number}

function CosmoCuboid(lowerleft::Vector{T}, upperright::Vector{T}) where {T<:Number}
    return CosmoCuboid{T}(lowerleft, upperright)
end

function CosmoCuboid(center::Vector{<:Number}, x::Number, y::Number, z::Number)
    halflengths = [x, y, z]
    halflengths *= 1 // 2
    return CosmoCuboid(center .- halflengths, center .+ halflengths)
end


"""
    CosmoRectangle

Alias for a 2D [`CosmoHyperrectangle`](@ref).
"""
const CosmoRectangle{T} = CosmoHyperrectangle{T,2} where {T<:Number}

function CosmoRectangle(lowerleft::Vector{T}, upperright::Vector{T}) where {T<:Number}
    return CosmoRectangle{T}(lowerleft, upperright)
end

function CosmoRectangle(center::Vector{<:Number}, x::Number, y::Number)
    halflengths = [x, y]
    halflengths *= 1 // 2
    return CosmoRectangle(center .- halflengths, center .+ halflengths)
end


"""
    CosmoHypercube(center::AbstractVector{<:Number}, radius::Number)

Return a cubic `CosmoCuboid` centered around `center`, with equal sidelengths.
"""
function CosmoHypercube(center::Vector{<:Number}, radius::Number)
    return CosmoHyperrectangle(center .- radius, center .+ radius)
end

"""
    CosmoCube(center::AbstractVector{<:Number}, radius::Number)

Return a cubic `CosmoCuboid` centered around `center`, with equal sidelengths.
"""
function CosmoCube(center::Vector{<:Number}, radius::Number)
    return CosmoCuboid(center .- radius, center .+ radius)
end

"""
    CosmoSquare(center::AbstractVector{<:Number}, radius::Number)

Return a square `CosmoRectangle` centered around `center`, with equal sidelengths.
"""
function CosmoSquare(center::Vector{<:Number}, radius::Number)
    return CosmoRectangle(center .- radius, center .+ radius)
end


for (name, N) in zip([:CosmoHyperrectangle, :CosmoCuboid, :CosmoRectangle], [:(length(lowerleft)), 3, 2])
    quote
        function $name(lowerleft::Vector{T1}, upperright::Vector{T2}) where {T1<:Number,T2<:Number}
            T1 === T2 && return CosmoHyperrectangle{T1,$N}(lowerleft, upperright)
            T = promote_type(T1, T2)
            return CosmoHyperrectangle{T,$N}(convert(Vector{T}, lowerleft), convert(Vector{T}, upperright))
        end
    end |> eval
end


geometry_enclosing_corners(rectangle::CosmoHyperrectangle) = rectangle.lowerleft, rectangle.upperright

function geometry_enclosing_center(rectangle::CosmoHyperrectangle)
    return 1 // 2 .* (rectangle.lowerleft .+ rectangle.upperright)
end

function mask_in(pos::AbstractMatrix{<:Number}, rectangle::CosmoCuboid{T}) where {T}
    @inbounds @views @. (
        (rectangle.lowerleft[1] ≤ pos[1, :] ≤ rectangle.upperright[1]) &
        (rectangle.lowerleft[2] ≤ pos[2, :] ≤ rectangle.upperright[2]) &
        (rectangle.lowerleft[3] ≤ pos[3, :] ≤ rectangle.upperright[3])
    )
end

function mask_in(pos::AbstractMatrix{<:Number}, rectangle::CosmoRectangle{T}) where {T}
    @inbounds @views @. (
        (rectangle.lowerleft[1] ≤ pos[1, :] ≤ rectangle.upperright[1]) &
        (rectangle.lowerleft[2] ≤ pos[2, :] ≤ rectangle.upperright[2])
    )
end

function mask_in(pos::AbstractMatrix{<:Number}, rectangle::CosmoHyperrectangle{T,N}) where {T,N}
    mask = @inbounds @views @. rectangle.lowerleft[1] ≤ pos[1, :] ≤ rectangle.upperright[1]
    for i in 2:N
        mask .&= @inbounds @views @. rectangle.lowerleft[i] ≤ pos[i, :] ≤ rectangle.upperright[i]
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

    function CosmoHypersphere{T,N}(center::Vector{T}, radius::T) where {T<:Number,N}
        @assert length(center) == N
        new{T,N}(center, radius)
    end
end

"""
    CosmoSphere

Alias for a 2D [`CosmoHypersphere`](@ref).
"""
const CosmoSphere{T} = CosmoHypersphere{T,3} where {T<:Number}

"""
    CosmoCircle

Alias for a 2D [`CosmoHypersphere`](@ref).
"""
const CosmoCircle{T} = CosmoHypersphere{T,2} where {T<:Number}

for (name, N) in zip([:CosmoHypersphere, :CosmoSphere, :CosmoCircle], [:(length(center)), 3, 2])
    quote
        function $name(center::Vector{T1}, radius::T2) where {T1<:Number,T2<:Number}
            T1 === T2 && return CosmoHypersphere{T1,$N}(center, radius)
            T = promote_type(T1, T2)
            return CosmoHypersphere{T,$N}(convert(Vector{T}, center), convert(T, radius))
        end
    end |> eval
end


function Base.:(==)(s1::CosmoHypersphere, s2::CosmoHypersphere)
    return s1.center == s2.center && s1.radius == s2.radius
end


function geometry_enclosing_corners(sphere::CosmoHypersphere)
    return sphere.center .- sphere.radius, sphere.center .+ sphere.radius
end

geometry_enclosing_center(sphere::CosmoHypersphere) = sphere.center

function mask_in(pos::AbstractMatrix{<:Number}, sphere::CosmoSphere)
    r² = sphere.radius^2
    @inbounds @views @. (pos[1, :] - sphere.center[1])^2 +
                        (pos[2, :] - sphere.center[2])^2 +
                        (pos[3, :] - sphere.center[3])^2 ≤ r²
end

function mask_in(pos::AbstractMatrix{<:Number}, circle::CosmoCircle)
    r² = circle.radius^2
    @inbounds @views @. (pos[1, :] - circle.center[1])^2 + (pos[2, :] - circle.center[2])^2 ≤ r²
end

function mask_in(pos::AbstractMatrix{<:Number}, sphere::CosmoHypersphere{T,N}) where {T,N}
    r² = sphere.radius^2
    center = sphere.center
    mask = BitArray(undef, size(pos, 2))
    @inbounds for i in eachindex(mask)
        s = (pos[1, i] - center[1])^2
        for j in 2:length(center)
            s += (pos[j, i] - center[j])^2
        end
        mask[i] = s ≤ r²
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

    function CosmoCylinder(startpos::Vector{T1}, endpos::Vector{T1}, radius::T2) where {T1<:Number,T2<:Number}
        T1 === T2 && return new{T1}(startpos, endpos, radius)
        T = promote_type(T1, T2)
        return new{T}(convert(Vector{T}, startpos), convert(Vector{T}, endpos), convert(T, radius))
    end
end

function Base.:(==)(c1::CosmoCylinder, c2::CosmoCylinder)
    return c1.startpos == c2.startpos && c1.endpos == c2.endpos && c1.radius == c2.radius
end


function geometry_enclosing_corners(cylinder::CosmoCylinder)
    min.(cylinder.startpos, cylinder.endpos) .- cylinder.radius,
    max.(cylinder.startpos, cylinder.endpos) .+ cylinder.radius
end

geometry_enclosing_center(cylinder::CosmoCylinder) = 1 // 2 .* (cylinder.startpos .+ cylinder.endpos)

function mask_in(pos::AbstractMatrix{<:Number}, cylinder::CosmoCylinder)
    # https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
    pos1 = @view pos[1, :]
    pos2 = @view pos[2, :]
    pos3 = @view pos[3, :]
    p1 = cylinder.startpos
    p2 = cylinder.endpos
    Δ = p2 .- p1
    norm2Δ = Δ[1]^2 + Δ[2]^2 + Δ[3]^2
    r² = cylinder.radius^2
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

    function CosmoStandingCylinder(
        center::Vector{T1},
        height::T2,
        radius::T3,
    ) where {T1<:Number,T2<:Number,T3<:Number}
        T1 === T2 === T3 && return new{T1}(center, height, radius)
        T = promote_type(T1, T2, T3)
        return new{T}(convert(Vector{T}, center), convert(T, height), convert(T, radius))
    end
end

function Base.:(==)(c1::CosmoStandingCylinder, c2::CosmoStandingCylinder)
    return c1.center == c2.center && c1.height == c2.height && c1.radius == c2.radius
end

function Base.:(==)(c1::CosmoCylinder, c2::CosmoStandingCylinder)
    return geometry_enclosing_center(c1) == geometry_enclosing_center(c2) &&
           abs.(c1.startpos[3] - c1.endpos[3]) == c2.height &&
           c1.radius == c2.radius &&
           @views c1.startpos[1:2] == c1.endpos[1:2] # c1 is actually standing
end
Base.:(==)(c1::CosmoStandingCylinder, c2::CosmoCylinder) = ==(c2, c1)


function geometry_enclosing_corners(cylinder::CosmoStandingCylinder)
    cylinder.center .- [cylinder.radius, cylinder.radius, 1 // 2 * cylinder.height],
    cylinder.center .+ [cylinder.radius, cylinder.radius, 1 // 2 * cylinder.height]
end

geometry_enclosing_center(cylinder::CosmoStandingCylinder) = cylinder.center

function mask_in(pos::AbstractMatrix{<:Number}, cylinder::CosmoStandingCylinder)
    r² = cylinder.radius^2
    halfheight = 1 // 2 * cylinder.height
    @inbounds @views @. ((pos[1, :] - cylinder.center[1])^2 + (pos[2, :] - cylinder.center[2])^2 ≤ r²) &
                        (cylinder.center[3] - halfheight ≤ pos[3, :] ≤ cylinder.center[3] + halfheight)
end
