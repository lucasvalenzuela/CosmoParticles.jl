@doc raw"""
    colnormell2(a::AbstractMatrix, q[, s])

Returns a new `Vector` with the columnwise squared elliptical norms of `a`.

In the 2D case: ``x^2 + \frac{y^2}{q^2}``.
In the 3D case: ``x^2 + \frac{y^2}{q^2} + \frac{z^2}{s^2}``.
"""
function colnormell2(a::AbstractMatrix{T}, q::S, _::Nothing=nothing) where {T,S}
    Tnew = float(promote_type(T, S))
    dst = Vector{Tnew}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, q)
end

function colnormell2(a::AbstractMatrix{Q}, q::S, _::Nothing=nothing) where {T,Q<:Quantity{T},S}
    Tnew = float(promote_type(T, S))
    u = unit(Q)^2
    Q² = Quantity{Tnew,dimension(u),typeof(u)}
    dst = Vector{Q²}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, q)
end

function colnormell2(a::AbstractMatrix{Union{Missing,Q}}, q::S, _::Nothing=nothing) where {T,Q<:Quantity{T},S}
    Tnew = float(promote_type(T, S))
    u = unit(Q)^2
    Q² = Quantity{Tnew,dimension(u),typeof(u)}
    dst = Vector{Union{Missing,Q²}}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, q)
end

function colnormell2(a::AbstractMatrix{T}, q::S, s::R) where {T,S,R}
    Tnew = float(promote_type(T, S, R))
    dst = Vector{Tnew}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, q, s)
end

function colnormell2(a::AbstractMatrix{Q}, q::S, s::R) where {T,Q<:Quantity{T},S,R}
    Tnew = float(promote_type(T, S, R))
    u = unit(Q)^2
    Q² = Quantity{Tnew,dimension(u),typeof(u)}
    dst = Vector{Q²}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, q, s)
end

function colnormell2(a::AbstractMatrix{Union{Missing,Q}}, q::S, s::R) where {T,Q<:Quantity{T},S,R}
    Tnew = float(promote_type(T, S, R))
    u = unit(Q)^2
    Q² = Quantity{Tnew,dimension(u),typeof(u)}
    dst = Vector{Union{Missing,Q²}}(undef, size(a, 2))
    return _colnorm2_unsafe!(dst, a, q, s)
end

function _colnorm2_unsafe!(dst::AbstractVector, a::AbstractMatrix, q)
    invq² = 1 / q^2
    @inbounds @simd for i in axes(a, 2)
        dst[i] = a[1, i]^2 + a[2, i]^2 * invq²
    end
    return dst
end

function _colnorm2_unsafe!(dst::AbstractVector, a::AbstractMatrix, q, s)
    invq², invs² = 1 / q^2, 1 / s^2
    @inbounds @simd for i in axes(a, 2)
        dst[i] = a[1, i]^2 + a[2, i]^2 * invq² + a[3, i]^2 * invs²
    end
    return dst
end



"""
    struct CosmoEllipse{T,U} <: AbstractCosmoGeometry where {T<:Number,U<:Real}
        a::T
        q::U
    end

Ellipse around the origin given by the semi-major axis `a` and the axis ratio `q = b/a` (where `b` is the semi-minor axis).
"""
struct CosmoEllipse{T,U} <: AbstractCosmoGeometry where {T<:Number,U<:Real}
    a::T
    q::U
end

"""
    CosmoEllipse(a::T; q::U=1, constant_area=false) where {T<:Number,U<:Real}

Returns an ellipse around the origin with semi-major axis `a` and the axis ratio `q`.

If `constant_area=true`, the semi-major axis of the ellipse is picked such that its area is equal to that of a
circle with radius `a`.

The axis ratio `q` has to be between 0 and 1.
"""
function CosmoEllipse(a::T; q::U=1, constant_area=false) where {T<:Number,U<:Real}
    @assert 0 ≤ q ≤ 1
    if constant_area
        anew = a / sqrt(q)
        CosmoEllipse{float(typeof(anew)),float(U)}(anew, q)
    else
        CosmoEllipse{float(T),float(U)}(a, q)
    end
end

function CosmoEllipse(a::T1, b::T2) where {T1<:Number,T2<:Number}
    @assert a ≥ b
    q = b / a
    Tnew = float(promote_type(T1, T2))
    return CosmoEllipse{float(Tnew),typeof(q)}(a, q)
end

function CosmoParticles.geometry_enclosing_corners(e::CosmoEllipse)
    a = e.a
    b = a * e.q
    return [-a, -b], [a, b]
end

CosmoParticles.geometry_enclosing_center(::CosmoEllipse{T}) where {T} = T[0, 0]

function CosmoParticles.mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, e::CosmoEllipse{T}) where {T}
    @assert size(pos, 1) == 2
    @assert size(pos, 2) == length(mask)
    r² = e.a^2
    if e.q == 1
        return @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 ≤ r²
    else
        invq² = 1 / e.q^2
        return @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 * invq² ≤ r²
    end
end



"""
    struct CosmoEllipsoid{T,U} <: AbstractCosmoGeometry where {T<:Number,U<:Real}
        a::T
        q::U
        s::U
    end

Ellipsoid around the origin given by the largest semi-major axis `a`, and the axis ratios `q = b/a` and `s = c/a`
(where `b` and `c` are the semi-minor axes with `b ≥ c`).
"""
struct CosmoEllipsoid{T,U} <: AbstractCosmoGeometry where {T<:Number,U<:Real}
    a::T
    q::U
    s::U
end


"""
    CosmoEllipsoid(a::Number; q::Real=1, s::Real=1, constant_volume=false)

Returns an ellipsoid around the origin with semi-major axis `a` and the axis ratios `q` and `s`.

If `constant_volume=true`, the semi-major axis of the ellipsoid is picked such that its volume is equal to that
of a sphere with radius `a`.

Both `q` and `s` have to be between 0 and 1, and `q ≥ s` has to hold.
"""
function CosmoEllipsoid(a; q=1, s=1, constant_volume=false)
    _create_cosmoellipsoid_qs(a, q, s, constant_volume)
end

function _create_cosmoellipsoid_qs(a::T, q::U1, s::U2, constant_volume) where {T<:Number,U1<:Real,U2<:Real}
    @assert q ≥ s
    @assert 0 ≤ q ≤ 1 && 0 ≤ s ≤ 1
    Unew = float(promote_type(U1, U2))
    if constant_volume
        anew = a / cbrt(q * s)
        CosmoEllipsoid{float(typeof(anew)),Unew}(anew, q, s)
    else
        CosmoEllipsoid{float(T),Unew}(a, q, s)
    end
end

"""
    CosmoEllipsoid(a::Number, b::Number, c::Number)

Returns an ellipsoid around the origin with semi axes `a`, `b`, and `c`, with `a ≥ b ≥ c`.

If `constant_volume=true`, the semi-major axis of the ellipsoid is picked such that its volume is equal to that
of a sphere with radius `a`.

If `c` is `nothing`, a [`CosmoEllipse`](@ref) is created instead.
"""
function CosmoEllipsoid(a::T1, b::T2, c::T3) where {T1<:Number,T2<:Number,T3<:Number}
    @assert a ≥ b ≥ c
    q, s = b / a, c / a
    Tnew = float(promote_type(T1, T2, T3))
    U = promote_type(typeof(q), typeof(s))
    return CosmoEllipsoid{float(Tnew),U}(a, q, s)
end

function CosmoParticles.geometry_enclosing_corners(e::CosmoEllipsoid)
    a = e.a
    b = a * e.q
    c = a * e.s
    return [-a, -b, -c], [a, b, c]
end

CosmoParticles.geometry_enclosing_center(::CosmoEllipsoid{T}) where {T} = T[0, 0, 0]

function CosmoParticles.mask_in!(
    mask::BitVector,
    pos::AbstractMatrix{<:Number},
    e::CosmoEllipsoid{T},
) where {T}
    @assert size(pos, 1) == 3
    @assert size(pos, 2) == length(mask)
    r² = e.a^2
    if e.q == 1 && e.s == 1
        return @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 + pos[3, :]^2 ≤ r²
    else
        invq² = 1 / e.q^2
        invs² = 1 / e.s^2
        return @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 * invq² + pos[3, :]^2 * invs² ≤ r²
    end
end

CosmoEllipsoid(a::T1, b::T2, c::Nothing) where {T1<:Number,T2<:Number} = CosmoEllipse(a, b)
function _create_cosmoellipsoid_qs(a::T, q::U1, s::Nothing, constant_volume) where {T<:Number,U1<:Real}
    CosmoEllipse(a; q, constant_area=constant_volume)
end



"""
    struct CosmoHomoeoid2D{T,U} <: AbstractCosmoGeometry where {T<:Number,U<:Real}
        a1::T
        a2::T
        q::U
    end

Homoeoid in two dimensions around the origin given by the inner and outer semi-major axes `a1` and `a2`, and the
axis ratio `q = b/a`.
"""
struct CosmoHomoeoid2D{T,U} <: AbstractCosmoGeometry where {T<:Number,U<:Real}
    a1::T
    a2::T
    q::U
end

"""
    CosmoHomoeoid2D(a1::Number, a2::Number; q::Real=1, constant_area=false)
    CosmoHomoeoid2D(a::Number, b::Number; a2::Number)

Returns a 2D homoeoid around the origin with semi-major axes `a1 ≤ a2` and the axis ratio `q`.

If `constant_area=true`, the semi-major axes of the ellipse are picked such that the enclosed area is equal to
that enclosed by two circles with radii `a1` and `a2`.

The axis ratio `q` has to be between 0 and 1.

For the second method, `a2` can also be smaller than `a1`. It is only important that `a` and `b` both either
correspond to the inner or outer ellipse.
"""
function CosmoHomoeoid2D(a1, a2_or_b; a2=nothing, q=1, constant_area=false)
    if isnothing(a2)
        _cosmo_homoeoid_2d(a1, a2_or_b, q, constant_area)
    else
        _cosmo_homoeoid_2d(a1, a2_or_b, a2)
    end
end

function _cosmo_homoeoid_2d(a1::T1, a2::T2, q::U, constant_area) where {T1<:Number,T2<:Number,U<:Real}
    @assert a1 < a2
    @assert 0 ≤ q ≤ 1
    if constant_area
        a1new = a1 / sqrt(q)
        a2new = a2 / sqrt(q)
        T = promote_type(float(typeof(a1new)), float(typeof(a2new)))
        CosmoHomoeoid2D{T,float(U)}(a1new, a2new, q)
    else
        T = promote_type(float(T1), float(T2))
        CosmoHomoeoid2D{T,float(U)}(a1, a2, q)
    end
end

function _cosmo_homoeoid_2d(a::T1, b::Number, a2::T2) where {T1<:Number,T2<:Number}
    @assert a ≥ b
    q = b / a
    a1, a2 = minmax(a, a2)
    Tnew = float(promote_type(T1, T2))
    return CosmoHomoeoid2D{Tnew,typeof(q)}(a1, a2, q)
end

function CosmoParticles.geometry_enclosing_corners(e::CosmoHomoeoid2D)
    a = e.a2
    b = a * e.q
    return [-a, -b], [a, b]
end

CosmoParticles.geometry_enclosing_center(::CosmoHomoeoid2D{T}) where {T} = T[0, 0]

function CosmoParticles.mask_in!(
    mask::BitVector,
    pos::AbstractMatrix{<:Number},
    e::CosmoHomoeoid2D{T},
) where {T}
    @assert size(pos, 1) == 2
    @assert size(pos, 2) == length(mask)
    r1² = e.a1^2
    r2² = e.a2^2
    if e.q == 1
        return @inbounds @views @. mask = r1² ≤ pos[1, :]^2 + pos[2, :]^2 ≤ r2²
    else
        invq² = 1 / e.q^2
        return @inbounds @views @. mask = r1² ≤ pos[1, :]^2 + pos[2, :]^2 * invq² ≤ r2²
    end
end



"""
    struct CosmoHomoeoid{T,U} <: AbstractCosmoGeometry where {T<:Number,U<:Real}
        a1::T
        a2::T
        q::U
        s::U
    end

Homoeoid in three dimensions around the origin given by the inner and outer semi-major axes `a1` and `a2`, and the
axis ratios `q = b/a` and `s = c/a`.
"""
struct CosmoHomoeoid{T,U} <: AbstractCosmoGeometry where {T<:Number,U<:Real}
    a1::T
    a2::T
    q::U
    s::U
end


"""
    CosmoHomoeoid(a1::Number, a2::Number; q::Real=1, s::Real=1, constant_volume=false)

Returns a 3D homoeoid around the origin with semi-major axes `a1 ≤ a2` and the axis ratios `q` and `s`.

If `constant_volume=true`, the semi-major axes of the ellipsoids are picked such that the enclosed volume is
equal to that enclosed by two spheres with radii `a1` and `a2`.

The axis ratios `q` and `s` have to be between 0 and 1.
"""
function CosmoHomoeoid(a1, a2; q=1, s=1, constant_volume=false)
    _create_cosmohomoeoid_qs(a1, a2, q, s, constant_volume)
end

function _create_cosmohomoeoid_qs(
    a1::T1,
    a2::T2,
    q::U1,
    s::U2,
    constant_volume,
) where {T1<:Number,T2<:Number,U1<:Real,U2<:Real}
    @assert q ≥ s
    @assert 0 ≤ q ≤ 1 && 0 ≤ s ≤ 1
    Unew = float(promote_type(U1, U2))
    if constant_volume
        factor = 1 / cbrt(q * s)
        a1new = a1 * factor
        a2new = a2 * factor
        Tnew = float(promote_type(typeof(a1new), typeof(a2new)))
        CosmoHomoeoid{Tnew,Unew}(a1new, a2new, q, s)
    else
        Tnew = float(promote_type(T1, T2))
        CosmoHomoeoid{float(Tnew),Unew}(a1, a2, q, s)
    end
end

"""
    CosmoHomoeoid(a::Number, b::Number, c::Number; a2::Number)

Returns a 3D homoeoid around the origin with semi axes `a ≥ b ≥ c` and another enclosing similar ellipsoid
specified by the semi-major axis `a2`.

The semi-major axis `a2` can be smaller or larger than `a`.

If `c` is `nothing`, a [`CosmoHomoeoid2D`](@ref) is created instead.
"""
function CosmoHomoeoid(a::T1, b::Number, c::Number; a2::T2) where {T1<:Number,T2<:Number}
    @assert a ≥ b ≥ c
    q, s = b / a, c / a
    a1, a2 = minmax(a, a2)
    Tnew = float(promote_type(T1, T2))
    U = promote_type(typeof(q), typeof(s))
    return CosmoHomoeoid{Tnew,U}(a1, a2, q, s)
end

function CosmoParticles.geometry_enclosing_corners(e::CosmoHomoeoid)
    a = e.a2
    b = a * e.q
    c = a * e.s
    return [-a, -b, -c], [a, b, c]
end

CosmoParticles.geometry_enclosing_center(::CosmoHomoeoid{T}) where {T} = T[0, 0, 0]

function CosmoParticles.mask_in!(
    mask::BitVector,
    pos::AbstractMatrix{<:Number},
    e::CosmoHomoeoid{T},
) where {T}
    @assert size(pos, 1) == 3
    @assert size(pos, 2) == length(mask)
    r1² = e.a1^2
    r2² = e.a2^2
    if e.q == 1 && e.s == 1
        return @inbounds @views @. mask = r1² ≤ pos[1, :]^2 + pos[2, :]^2 + pos[3, :]^2 ≤ r2²
    else
        invq² = 1 / e.q^2
        invs² = 1 / e.s^2
        return @inbounds @views @. mask = r1² ≤ pos[1, :]^2 + pos[2, :]^2 * invq² + pos[3, :]^2 * invs² ≤ r2²
    end
end

CosmoHomoeoid(a::T1, b::T2, c::Nothing; a2) where {T1<:Number,T2<:Number} = CosmoHomoeoid2D(a, b; a2)
function _create_cosmohomoeoid_qs(
    a1::T1,
    a2::T2,
    q::U1,
    s::Nothing,
    constant_volume,
) where {T1<:Number,T2<:Number,U1<:Real}
    CosmoHomoeoid2D(a1, a2; q, constant_area=constant_volume)
end


"""
    triaxiality(q::Real, s::Real)
    triaxiality(e::Union{CosmoEllipsoid,CosmoHomoeoid})

Triaxiality from the axis ratios `q` and `s`.

-   ``0 < T < 1/3``: oblate
- ``1/3 < T < 2/3``: triaxial
- ``2/3 < T < 1``:   prolate
"""
@inline triaxiality(q::Real, s::Real) = (1 - q^2) / (1 - s^2)
@inline triaxiality(e::Union{CosmoEllipsoid,CosmoHomoeoid}) = triaxiality(e.q, e.s)

@doc raw"""
    ellipticity(q::Real)
    ellipticity(e::Union{CosmoEllipse,CosmoHomoeoid2D})

Ellipticity from axis ratio `q`: ``ϵ = 1 - q``
"""
@inline ellipticity(q::Real) = 1 - q
@inline ellipticity(e::Union{CosmoEllipse,CosmoHomoeoid2D}) = ellipticity(e.q)

@doc raw"""
    eccentricity(q::Real)
    eccentricity(e::Union{CosmoEllipse,CosmoHomoeoid2D})

Eccentricity from axis ratio `q`: ``e = \sqrt{1 - q^2}``
"""
@inline eccentricity(q::Real) = sqrt(1 - q^2)
@inline eccentricity(e::Union{CosmoEllipse,CosmoHomoeoid2D}) = eccentricity(e.q)
