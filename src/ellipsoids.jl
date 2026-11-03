@doc raw"""
    colnormell2(a::AbstractMatrix, q[, s])

Returns a new `Vector` with the columnwise squared elliptical norms of `a`.

In the 2D case: ``x^2 + \frac{y^2}{q^2}``.
In the 3D case: ``x^2 + \frac{y^2}{q^2} + \frac{z^2}{s^2}``.

See also [`colnormell`](@ref).

# Examples

```jldoctest
julia> a = [3.0 0.0 0.0; 0.0 2.0 0.5];

julia> colnormell2(a, 0.5)
3-element Vector{Float64}:
  9.0
 16.0
  1.0
```
"""
function colnormell2(a::AbstractMatrix{T}, q, ::Nothing=nothing) where {T}
    # the floating-point type of the elements of `a` promoted with the axis ratios
    Tnew = Base.promote_op((x, q) -> float(x^2 + x^2 * q), T, typeof(q))
    dst = Vector{Tnew}(undef, size(a, 2))
    return _colnormell2_unsafe!(dst, a, q)
end

function colnormell2(a::AbstractMatrix{T}, q, s) where {T}
    # the floating-point type of the elements of `a` promoted with the axis ratios
    Tnew = Base.promote_op((x, q, s) -> float(x^2 + x^2 * q + x^2 * s), T, typeof(q), typeof(s))
    dst = Vector{Tnew}(undef, size(a, 2))
    return _colnormell2_unsafe!(dst, a, q, s)
end

# The elliptical-norm kernels index the rows of `a` from 1: two rows for `q`, three for `q` and `s`.
function _check_ellnorm_rows(a::AbstractMatrix, nrows)
    Base.require_one_based_indexing(a)
    size(a, 1) == nrows || throw(
        DimensionMismatch(
            "expected a matrix with $nrows rows for the given axis ratios, got $(size(a, 1)) rows",
        ),
    )
    return nothing
end

function _colnormell2_unsafe!(dst::AbstractVector, a::AbstractMatrix, q)
    _check_ellnorm_rows(a, 2)
    invq² = 1 / q^2
    @inbounds @simd for i in axes(a, 2)
        dst[i] = a[1, i]^2 + a[2, i]^2 * invq²
    end
    return dst
end

function _colnormell2_unsafe!(dst::AbstractVector, a::AbstractMatrix, q, s)
    _check_ellnorm_rows(a, 3)
    invq², invs² = 1 / q^2, 1 / s^2
    @inbounds @simd for i in axes(a, 2)
        dst[i] = a[1, i]^2 + a[2, i]^2 * invq² + a[3, i]^2 * invs²
    end
    return dst
end



@doc raw"""
    colnormell(a::AbstractMatrix, q[, s])

Returns a new `Vector` with the columnwise elliptical norms of `a`.

In the 2D case: ``\sqrt{x^2 + \frac{y^2}{q^2}}``.
In the 3D case: ``\sqrt{x^2 + \frac{y^2}{q^2} + \frac{z^2}{s^2}}``.

See also [`colnormell2`](@ref).

# Examples

```jldoctest
julia> a = [3.0 0.0 0.0; 0.0 2.0 0.5];

julia> colnormell(a, 0.5)
3-element Vector{Float64}:
 3.0
 4.0
 1.0
```
"""
function colnormell(a::AbstractMatrix{T}, q, ::Nothing=nothing) where {T}
    # the floating-point type of the elements of `a` promoted with the axis ratios
    Tnew = Base.promote_op((x, q) -> sqrt(x^2 + x^2 * q), T, typeof(q))
    dst = Vector{Tnew}(undef, size(a, 2))
    return _colnormell_unsafe!(dst, a, q)
end

function colnormell(a::AbstractMatrix{T}, q, s) where {T}
    # the floating-point type of the elements of `a` promoted with the axis ratios
    Tnew = Base.promote_op((x, q, s) -> sqrt(x^2 + x^2 * q + x^2 * s), T, typeof(q), typeof(s))
    dst = Vector{Tnew}(undef, size(a, 2))
    return _colnormell_unsafe!(dst, a, q, s)
end

function _colnormell_unsafe!(dst::AbstractVector, a::AbstractMatrix, q)
    _check_ellnorm_rows(a, 2)
    invq² = 1 / q^2
    @inbounds @simd for i in axes(a, 2)
        dst[i] = sqrt(a[1, i]^2 + a[2, i]^2 * invq²)
    end
    return dst
end

function _colnormell_unsafe!(dst::AbstractVector, a::AbstractMatrix, q, s)
    _check_ellnorm_rows(a, 3)
    invq², invs² = 1 / q^2, 1 / s^2
    @inbounds @simd for i in axes(a, 2)
        dst[i] = sqrt(a[1, i]^2 + a[2, i]^2 * invq² + a[3, i]^2 * invs²)
    end
    return dst
end



function _check_axis_ratios(q)
    0 < q ≤ 1 || throw(ArgumentError("the axis ratio q must satisfy 0 < q ≤ 1, got q = $q"))
    return nothing
end

function _check_axis_ratios(q, s)
    _check_axis_ratios(q)
    0 < s ≤ q || throw(ArgumentError("the axis ratio s must satisfy 0 < s ≤ q, got q = $q and s = $s"))
    return nothing
end

function _check_semimajor_axes(a1, a2)
    a1 ≤ a2 || throw(
        ArgumentError(
            "the inner semi-major axis a1 must not exceed the outer one a2, got a1 = $a1 and a2 = $a2",
        ),
    )
    return nothing
end

"""
    struct CosmoEllipse{T<:Number,U<:Real} <: AbstractCosmoGeometry
        a::T
        q::U
    end

Ellipse around the origin given by the semi-major axis `a` and the axis ratio `q = b/a` (where `b` is the semi-minor axis).

The fields store the semi-major axis and the axis ratio; construct it with [`CosmoEllipse(a, b)`](@ref CosmoEllipse)
from the semi-axes or [`CosmoEllipse(a; q)`](@ref CosmoEllipse) from the axis ratio.

# Examples

```jldoctest
julia> e = CosmoEllipse(2.0, 1.0)
CosmoEllipse{Float64, Float64}(2.0, 0.5)

julia> CosmoParticles.mask_in([1.0 0.0 2.0; 0.0 1.5 0.5], e)
3-element BitVector:
 1
 0
 0
```
"""
struct CosmoEllipse{T<:Number,U<:Real} <: AbstractCosmoGeometry
    a::T
    q::U

    function CosmoEllipse{T,U}(a, q) where {T,U}
        e = new{T,U}(a, q)
        _check_axis_ratios(e.q)
        return e
    end
end

"""
    CosmoEllipse(a::T; q::U=1, constant_area=false) where {T<:Number,U<:Real}

Returns an ellipse around the origin with semi-major axis `a` and the axis ratio `q`.

If `constant_area=true`, the semi-major axis of the ellipse is picked such that its area is equal to that of a
circle with radius `a`.

The axis ratio `q` has to satisfy `0 < q ≤ 1`.

# Examples

```jldoctest
julia> CosmoEllipse(2.0; q=0.5)
CosmoEllipse{Float64, Float64}(2.0, 0.5)
```
"""
function CosmoEllipse(a::T; q::U=1, constant_area=false) where {T<:Number,U<:Real}
    if constant_area
        _check_axis_ratios(q) # before `sqrt`, which throws for negative q
        anew = a / sqrt(q)
        CosmoEllipse{float(typeof(anew)),float(U)}(anew, q)
    else
        CosmoEllipse{float(T),float(U)}(a, q)
    end
end

"""
    CosmoEllipse(a::Number, b::Number)

Returns an ellipse around the origin with semi-major axis `a` and semi-minor axis `b`.

The semi-minor axis has to be smaller than the semi-major axis.

# Examples

```jldoctest
julia> CosmoEllipse(2.0, 1.0)
CosmoEllipse{Float64, Float64}(2.0, 0.5)
```
"""
function CosmoEllipse(a::T1, b::T2) where {T1<:Number,T2<:Number}
    q = uconvert(NoUnits, b / a)
    Tnew = float(promote_type(T1, T2))
    return CosmoEllipse{float(Tnew),typeof(q)}(a, q)
end

function geometry_enclosing_corners(e::CosmoEllipse)
    a = e.a
    b = a * e.q
    return [-a, -b], [a, b]
end

geometry_enclosing_center(::CosmoEllipse{T}) where {T} = zeros(T, 2)

function mask_in!(mask::BitVector, pos::AbstractMatrix{<:Number}, e::CosmoEllipse{T}) where {T}
    _check_mask_dims(mask, pos, 2)
    Base.require_one_based_indexing(pos)
    r² = e.a^2
    if e.q == 1
        return @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 ≤ r²
    else
        invq² = 1 / e.q^2
        return @inbounds @views @. mask = pos[1, :]^2 + pos[2, :]^2 * invq² ≤ r²
    end
end



"""
    struct CosmoEllipsoid{T<:Number,U<:Real} <: AbstractCosmoGeometry
        a::T
        q::U
        s::U
    end

Ellipsoid around the origin given by the largest semi-major axis `a`, and the axis ratios `q = b/a` and `s = c/a`
(where `b` and `c` are the semi-minor axes with `b ≥ c`).

The fields store the semi-major axis and the axis ratios; construct it with
[`CosmoEllipsoid(a, b, c)`](@ref CosmoEllipsoid) from the semi-axes or
[`CosmoEllipsoid(a; q, s)`](@ref CosmoEllipsoid) from the axis ratios.

# Examples

```jldoctest
julia> e = CosmoEllipsoid(4.0, 2.0, 1.0)
CosmoEllipsoid{Float64, Float64}(4.0, 0.5, 0.25)

julia> CosmoParticles.mask_in([3.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 0.5], e)
3-element BitVector:
 1
 0
 1
```
"""
struct CosmoEllipsoid{T<:Number,U<:Real} <: AbstractCosmoGeometry
    a::T
    q::U
    s::U

    function CosmoEllipsoid{T,U}(a, q, s) where {T,U}
        e = new{T,U}(a, q, s)
        _check_axis_ratios(e.q, e.s)
        return e
    end
end


"""
    CosmoEllipsoid(a::Number; q::Real=1, s::Real=q, constant_volume=false)

Returns an ellipsoid around the origin with semi-major axis `a` and the axis ratios `q` and `s`.

If `constant_volume=true`, the semi-major axis of the ellipsoid is picked such that its volume is equal to that
of a sphere with radius `a`.

The axis ratios have to satisfy `0 < s ≤ q ≤ 1`. By default, `s = q`: passing only `q` gives a prolate shape
and passing only `s` an oblate one.

# Examples

```jldoctest
julia> CosmoEllipsoid(4.0; q=0.5, s=0.25)
CosmoEllipsoid{Float64, Float64}(4.0, 0.5, 0.25)
```
"""
function CosmoEllipsoid(a; q=1, s=q, constant_volume=false)
    _create_cosmoellipsoid_qs(a, q, s, constant_volume)
end

function _create_cosmoellipsoid_qs(a::T, q::U1, s::U2, constant_volume) where {T<:Number,U1<:Real,U2<:Real}
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

If `c` is `nothing`, a [`CosmoEllipse`](@ref) is created instead.

# Examples

```jldoctest
julia> CosmoEllipsoid(4.0, 2.0, 1.0)
CosmoEllipsoid{Float64, Float64}(4.0, 0.5, 0.25)
```
"""
function CosmoEllipsoid(a::T1, b::T2, c::T3) where {T1<:Number,T2<:Number,T3<:Number}
    q, s = uconvert(NoUnits, b / a), uconvert(NoUnits, c / a)
    Tnew = float(promote_type(T1, T2, T3))
    U = promote_type(typeof(q), typeof(s))
    return CosmoEllipsoid{float(Tnew),U}(a, q, s)
end

function geometry_enclosing_corners(e::CosmoEllipsoid)
    a = e.a
    b = a * e.q
    c = a * e.s
    return [-a, -b, -c], [a, b, c]
end

geometry_enclosing_center(::CosmoEllipsoid{T}) where {T} = zeros(T, 3)

function mask_in!(
    mask::BitVector,
    pos::AbstractMatrix{<:Number},
    e::CosmoEllipsoid{T},
) where {T}
    _check_mask_dims(mask, pos, 3)
    Base.require_one_based_indexing(pos)
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
    struct CosmoHomoeoid2D{T<:Number,U<:Real} <: AbstractCosmoGeometry
        a1::T
        a2::T
        q::U
    end

Homoeoid in two dimensions around the origin given by the inner and outer semi-major axes `a1` and `a2`, and the
axis ratio `q = b/a`.

The fields store the semi-major axes and the axis ratio; construct it with
[`CosmoHomoeoid2D(a, b; a2)`](@ref CosmoHomoeoid2D) from the semi-axes of one ellipse and the other semi-major axis,
or [`CosmoHomoeoid2D(a1, a2; q)`](@ref CosmoHomoeoid2D) from the axis ratio.

# Examples

```jldoctest
julia> h = CosmoHomoeoid2D(1.0, 2.0; q=0.5)
CosmoHomoeoid2D{Float64, Float64}(1.0, 2.0, 0.5)

julia> CosmoParticles.mask_in([0.5 1.5 3.0; 0.0 0.0 0.0], h)
3-element BitVector:
 0
 1
 0
```
"""
struct CosmoHomoeoid2D{T<:Number,U<:Real} <: AbstractCosmoGeometry
    a1::T
    a2::T
    q::U

    function CosmoHomoeoid2D{T,U}(a1, a2, q) where {T,U}
        e = new{T,U}(a1, a2, q)
        _check_semimajor_axes(e.a1, e.a2)
        _check_axis_ratios(e.q)
        return e
    end
end

"""
    CosmoHomoeoid2D(a1::Number, a2::Number; q::Real=1, constant_area=false)
    CosmoHomoeoid2D(a::Number, b::Number; a2::Number)

Returns a 2D homoeoid around the origin with semi-major axes `a1 ≤ a2` and the axis ratio `q`.

If `constant_area=true`, the semi-major axes of the ellipse are picked such that the enclosed area is equal to
that enclosed by two circles with radii `a1` and `a2`.

The axis ratio `q` has to satisfy `0 < q ≤ 1`.

For the second method, `a2` can also be smaller than `a1`. It is only important that `a` and `b` both either
correspond to the inner or outer ellipse.

# Examples

```jldoctest
julia> CosmoHomoeoid2D(1.0, 2.0; q=0.5)
CosmoHomoeoid2D{Float64, Float64}(1.0, 2.0, 0.5)

julia> CosmoHomoeoid2D(2.0, 1.0; a2=4.0)
CosmoHomoeoid2D{Float64, Float64}(2.0, 4.0, 0.5)
```
"""
function CosmoHomoeoid2D(a1, a2_or_b; a2=nothing, q=1, constant_area=false)
    if isnothing(a2)
        _cosmo_homoeoid_2d(a1, a2_or_b, q, constant_area)
    else
        _cosmo_homoeoid_2d(a1, a2_or_b, a2)
    end
end

function _cosmo_homoeoid_2d(a1::T1, a2::T2, q::U, constant_area) where {T1<:Number,T2<:Number,U<:Real}
    if constant_area
        _check_axis_ratios(q) # before `sqrt`, which throws for negative q
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
    q = uconvert(NoUnits, b / a)
    a1, a2 = minmax(a, a2)
    Tnew = float(promote_type(T1, T2))
    return CosmoHomoeoid2D{Tnew,typeof(q)}(a1, a2, q)
end

function geometry_enclosing_corners(e::CosmoHomoeoid2D)
    a = e.a2
    b = a * e.q
    return [-a, -b], [a, b]
end

geometry_enclosing_center(::CosmoHomoeoid2D{T}) where {T} = zeros(T, 2)

function mask_in!(
    mask::BitVector,
    pos::AbstractMatrix{<:Number},
    e::CosmoHomoeoid2D{T},
) where {T}
    _check_mask_dims(mask, pos, 2)
    Base.require_one_based_indexing(pos)
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
    struct CosmoHomoeoid{T<:Number,U<:Real} <: AbstractCosmoGeometry
        a1::T
        a2::T
        q::U
        s::U
    end

Homoeoid in three dimensions around the origin given by the inner and outer semi-major axes `a1` and `a2`, and the
axis ratios `q = b/a` and `s = c/a`.

The fields store the semi-major axes and the axis ratios; construct it with
[`CosmoHomoeoid(a, b, c; a2)`](@ref CosmoHomoeoid) from the semi-axes of one ellipsoid and the other semi-major
axis, or [`CosmoHomoeoid(a1, a2; q, s)`](@ref CosmoHomoeoid) from the axis ratios.

# Examples

```jldoctest
julia> h = CosmoHomoeoid(1.0, 2.0; q=0.5, s=0.25)
CosmoHomoeoid{Float64, Float64}(1.0, 2.0, 0.5, 0.25)

julia> CosmoParticles.mask_in([0.5 1.5 3.0; 0.0 0.0 0.0; 0.0 0.0 0.0], h)
3-element BitVector:
 0
 1
 0
```
"""
struct CosmoHomoeoid{T<:Number,U<:Real} <: AbstractCosmoGeometry
    a1::T
    a2::T
    q::U
    s::U

    function CosmoHomoeoid{T,U}(a1, a2, q, s) where {T,U}
        e = new{T,U}(a1, a2, q, s)
        _check_semimajor_axes(e.a1, e.a2)
        _check_axis_ratios(e.q, e.s)
        return e
    end
end


"""
    CosmoHomoeoid(a1::Number, a2::Number; q::Real=1, s::Real=q, constant_volume=false)

Returns a 3D homoeoid around the origin with semi-major axes `a1 ≤ a2` and the axis ratios `q` and `s`.

If `constant_volume=true`, the semi-major axes of the ellipsoids are picked such that the enclosed volume is
equal to that enclosed by two spheres with radii `a1` and `a2`.

The axis ratios have to satisfy `0 < s ≤ q ≤ 1`. By default, `s = q`: passing only `q` gives a prolate shape
and passing only `s` an oblate one.

# Examples

```jldoctest
julia> CosmoHomoeoid(1.0, 2.0; q=0.5, s=0.25)
CosmoHomoeoid{Float64, Float64}(1.0, 2.0, 0.5, 0.25)
```
"""
function CosmoHomoeoid(a1, a2; q=1, s=q, constant_volume=false)
    _create_cosmohomoeoid_qs(a1, a2, q, s, constant_volume)
end

function _create_cosmohomoeoid_qs(
    a1::T1,
    a2::T2,
    q::U1,
    s::U2,
    constant_volume,
) where {T1<:Number,T2<:Number,U1<:Real,U2<:Real}
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

# Examples

```jldoctest
julia> CosmoHomoeoid(2.0, 1.0, 0.5; a2=4.0)
CosmoHomoeoid{Float64, Float64}(2.0, 4.0, 0.5, 0.25)
```
"""
function CosmoHomoeoid(a::T1, b::Number, c::Number; a2::T2) where {T1<:Number,T2<:Number}
    q, s = uconvert(NoUnits, b / a), uconvert(NoUnits, c / a)
    a1, a2 = minmax(a, a2)
    Tnew = float(promote_type(T1, T2))
    U = promote_type(typeof(q), typeof(s))
    return CosmoHomoeoid{Tnew,U}(a1, a2, q, s)
end

function geometry_enclosing_corners(e::CosmoHomoeoid)
    a = e.a2
    b = a * e.q
    c = a * e.s
    return [-a, -b, -c], [a, b, c]
end

geometry_enclosing_center(::CosmoHomoeoid{T}) where {T} = zeros(T, 3)

function mask_in!(
    mask::BitVector,
    pos::AbstractMatrix{<:Number},
    e::CosmoHomoeoid{T},
) where {T}
    _check_mask_dims(mask, pos, 3)
    Base.require_one_based_indexing(pos)
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

# Examples

```jldoctest
julia> triaxiality(0.5, 0.25)
0.8

julia> triaxiality(CosmoEllipsoid(4.0, 2.0, 1.0))
0.8
```
"""
@inline triaxiality(q::Real, s::Real) = (1 - q^2) / (1 - s^2)
@inline triaxiality(e::Union{CosmoEllipsoid,CosmoHomoeoid}) = triaxiality(e.q, e.s)

@doc raw"""
    ellipticity(q::Real)
    ellipticity(e::Union{CosmoEllipse,CosmoHomoeoid2D})

Ellipticity from axis ratio `q`: ``ϵ = 1 - q``

# Examples

```jldoctest
julia> ellipticity(0.5)
0.5

julia> ellipticity(CosmoEllipse(2.0, 1.0))
0.5
```
"""
@inline ellipticity(q::Real) = 1 - q
@inline ellipticity(e::Union{CosmoEllipse,CosmoHomoeoid2D}) = ellipticity(e.q)

@doc raw"""
    eccentricity(q::Real)
    eccentricity(e::Union{CosmoEllipse,CosmoHomoeoid2D})

Eccentricity from axis ratio `q`: ``e = \sqrt{1 - q^2}``

# Examples

```jldoctest
julia> eccentricity(0.5)
0.8660254037844386

julia> eccentricity(CosmoEllipse(2.0, 1.0))
0.8660254037844386
```
"""
@inline eccentricity(q::Real) = sqrt(1 - q^2)
@inline eccentricity(e::Union{CosmoEllipse,CosmoHomoeoid2D}) = eccentricity(e.q)

"""
    axis_lengths(e::Union{CosmoEllipse,CosmoEllipsoid})
    axis_lengths(e::Union{CosmoHomoeoid2D,CosmoHomoeoid}, region=:outer)

Returns the axis lengths of the ellipse or ellipsoid as a tuple.

For homoeoids, `region` selects the outer (`:outer`) or inner (`:inner`) boundary.

# Examples

```jldoctest
julia> axis_lengths(CosmoEllipse(4.0; q=0.5))
(4.0, 2.0)

julia> axis_lengths(CosmoEllipsoid(4.0; q=0.5, s=0.25))
(4.0, 2.0, 1.0)

julia> axis_lengths(CosmoHomoeoid(2.0, 4.0; q=0.5, s=0.25), :inner)
(2.0, 1.0, 0.5)
```
"""
axis_lengths(e::CosmoEllipse) = (e.a, e.a * e.q)
axis_lengths(e::CosmoEllipsoid) = (e.a, e.a * e.q, e.a * e.s)
axis_lengths(e::CosmoHomoeoid2D, region=:outer) = (a = _homoeoid_axis(e, region); (a, a * e.q))
axis_lengths(e::CosmoHomoeoid, region=:outer) = (a = _homoeoid_axis(e, region); (a, a * e.q, a * e.s))

function _homoeoid_axis(e, region)
    region === :outer && return e.a2
    region === :inner && return e.a1
    throw(ArgumentError("the region must be :inner or :outer, got $(repr(region))"))
end

for G in (CosmoEllipse, CosmoEllipsoid, CosmoHomoeoid2D, CosmoHomoeoid)
    @eval Base.:(==)(e1::$G, e2::$G) = _fields_equal(e1, e2)
    @eval Base.hash(e::$G, h::UInt) = _fields_hash(e, h)
end
