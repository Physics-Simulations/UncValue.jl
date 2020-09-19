module UncValue

import Base: +, -, *, /, \, ÷, ^, ==, !=, <, >, <=, >=
import Printf: @printf

struct Value{T <: Real}
    x::T
    σ::T

    Value{T}(x, σ) where {T <: Real} = σ < 0 ? error("Negative uncertainty has no meaning.") : new{T}(x, σ)
end
function Value(x::T, σ::U) where {T <: Real, U <: Real}
    xT, σT = promote(x, σ)
    @warn "Value was initialize with different types, promoting $(typeof(x)) & $(typeof(σ)) to $(typeof(xT))."
    Value{typeof(xT)}(xT, σT)
end
Value(x::T, σ::T) where {T <: Real} = Value{T}(x, σ)
Value(x::T) where {T <: Real} = Value{T}(x, zero(x))

export Value, val, unc, set_unc

## MATHEMATICAL OPERATORS

+(a::Value, b::T) where {T <: Real} = Value(a.x + b, a.σ)
+(a::T, b::Value) where {T <: Real} = Value(a + b.x, b.σ)
+(a::Value, b::Value) = Value(a.x + b.x, hypot(a.σ, b.σ))

-(a::Value) = Value(-a.x, a.σ)
-(a::Value, b::T) where {T <: Real} = Value(a.x - b, a.σ)
-(a::T, b::Value) where {T <: Real} = Value(a - b.x, b.σ)
-(a::Value, b::Value) = Value(a.x - b.x, hypot(a.σ, b.σ))

*(a::Value, b::T) where {T <: Real} = Value(a.x * b, a.σ * abs(b))
*(a::T, b::Value) where {T <: Real} = Value(a * b.x, b.σ * abs(a))
*(a::Value, b::Value) = Value(a.x * b.x, hypot(a.σ * b.x, a.x * b.σ))

/(a::Value, b::T) where {T <: Real} = Value(a.x / b, a.σ / abs(b))
/(a::T, b::Value) where {T <: Real} = Value(a / b.x, abs(a) * b.σ / b.x^2)
/(a::Value, b::Value) = Value(a.x / b.x, hypot(a.σ / b.x, b.σ * a.x / b.x^2))

\(a::Value, b::T) where {T <: Real} = Value(a.x \ b, a.x^2 \ b * a.σ)
\(a::T, b::Value) where {T <: Real} = Value(a \ b.x, a \ b.σ)
\(a::Value, b::Value) = Value(a.x \ b.x, hypot( a.x^2 \ b.x * a.σ, a.x \ b.σ))

÷(a::Value, b::T) where {T <: Real} = Value(a.x ÷ b, a.σ ÷ b)
÷(a::T, b::Value) where {T <: Real} = Value(a ÷ b.x, b.σ * a ÷ b.x^2)
÷(a::Value, b::Value) = Value(a.x ÷ b.x, round(typeof(a.x), hypot(a.σ ÷ b.x, b.σ * a.x ÷ b.x^2)))

^(a::Value, b::T) where {T <: Real} = Value(a.x ^ b, abs(b * a.x^(b-1)) * a.σ)
^(a::T, b::Value) where {T <: Real} = Value(a ^ b.x, abs(a^b.x * log(abs(a))) * b.σ)
^(a::Value, b::Value) where {T <: Real} =
    Value(a.x ^ b.x, hypot(b.x * a.x^(b.x-1) * a.σ, a.x^b.x * log(abs(a.x)) * b.σ))

==(a::Value, b::T) where {T <: Real} = a.x == b
==(a::T, b::Value) where {T <: Real} = a == b.x
==(a::Value, b::Value) = a.x == b.x

!=(a::Value, b::T) where {T <: Real} = a.x != b
!=(a::T, b::Value) where {T <: Real} = a != b.x
!=(a::Value, b::Value) = a.x != b.x

>(a::Value, b::T) where {T <: Real} = a.x > b
>(a::T, b::Value) where {T <: Real} = a > b.x
>(a::Value, b::Value) = a.x > b.x

<(a::Value, b::T) where {T <: Real} = a.x < b
<(a::T, b::Value) where {T <: Real} = a < b.x
<(a::Value, b::Value) = a.x < b.x

>=(a::Value, b::T) where {T <: Real} = a.x >= b
>=(a::T, b::Value) where {T <: Real} = a >= b.x
>=(a::Value, b::Value) = a.x >= b.x

<=(a::Value, b::T) where {T <: Real} = a.x <= b
<=(a::T, b::Value) where {T <: Real} = a <= b.x
<=(a::Value, b::Value) = a.x <= b.x

## MATHEMATICAL FUNCTIONS
Base.abs(a::Value) = Value(abs(a.x), a.σ)
Base.abs2(a::Value) = Value(abs2(a.x), 2 * abs(a.x) * a.σ)

Base.sin(a::Value) = Value(sin(a.x), abs(cos(a.x)) * a.σ)
Base.sind(a::Value) = Value(sind(a.x), abs(cosd(a.x)) * a.σ)
Base.sinh(a::Value) = Value(sinh(a.x), cosh(a.x) * a.σ)
Base.Math.sinpi(a::Value) = Value(sinpi(a.x), abs(π * cospi(a.x)) * a.σ)

Base.cos(a::Value) = Value(cos(a.x), abs(sin(a.x)) * a.σ)
Base.cosd(a::Value) = Value(cosd(a.x), abs(sind(a.x)) * a.σ)
Base.cosh(a::Value) = Value(cosh(a.x), abs(sinh(a.x)) * a.σ)
Base.Math.cospi(a::Value) = Value(cospi(a.x), abs(π * sinpi(a.x)) * a.σ)

Base.tan(a::Value) = Value(tan(a.x),  a.σ * sec(a.x)^2)
Base.tand(a::Value) = Value(tand(a.x),  a.σ * secd(a.x)^2)
Base.tanh(a::Value) = Value(tanh(a.x),  a.σ * sech(a.x)^2)

Base.asin(a::Value) = Value(asin(a.x),  a.σ / sqrt(1-a.x^2))
Base.Math.asind(a::Value) = Value(asind(a.x),  a.σ / sqrt(1-a.x^2))
Base.asinh(a::Value) = Value(asinh(a.x),  a.σ / sqrt(1 + a.x^2))

Base.acos(a::Value) = Value(acos(a.x),  a.σ / sqrt(1-a.x^2))
Base.Math.acosd(a::Value) = Value(acosd(a.x),  a.σ / sqrt(1-a.x^2))
Base.acosh(a::Value) = Value(acosh(a.x), a.σ / sqrt(a.x^2-1))

Base.atan(a::Value) = Value(atan(a.x), a.σ / (1 + a.x^2))
Base.atan(a::Value, b::T) where {T <: Real} = Value(atan(a.x, b), abs(b) * a.σ / (a.x^2 + b^2))
Base.atan(a::T, b::Value) where {T <: Real} = Value(atan(a, b.x), abs(a) * b.σ / (a^2 + b.x^2))
Base.atan(a::Value, b::Value) =
    Value(atan(a.x, b.x), hypot(b.x * a.σ, a.x * b.σ) / (a.x^2 + b.x^2))
Base.Math.atand(a::Value) = Value(atand(a.x),  a.σ / (1 + a.x^2))
Base.Math.atand(a::Value, b::T) where {T <: Real} = Value(atand(a.x, b), abs(b) * a.σ / (a.x^2 + b^2))
Base.Math.atand(a::T, b::Value) where {T <: Real} = Value(atand(a, b.x), abs(a) * b.σ / (a^2 + b.x^2))
Base.Math.atand(a::Value, b::Value) =
    Value(atand(a.x, b.x), hypot(b.x * a.σ, a.x * b.σ) / (a.x^2 + b.x^2))
Base.atanh(a::Value) = Value(atanh(a.x),  a.σ / (1 - a.x^2))

Base.Math.sec(a::Value) = Value(sec(a.x), abs(sec(a.x) * tan(a.x)) * a.σ)
Base.Math.secd(a::Value) = Value(secd(a.x), abs(secd(a.x) * tand(a.x)) * a.σ)
Base.Math.sech(a::Value) = Value(sech(a.x), abs(tanh(a.x) * sech(a.x)) * a.σ)

Base.Math.csc(a::Value) = Value(csc(a.x), abs(cot(a.x) * csc(a.x)) * a.σ)
Base.Math.cscd(a::Value) = Value(cscd(a.x), abs(cotd(a.x) * cscd(a.x)) * a.σ)
Base.Math.csch(a::Value) = Value(csch(a.x), abs(coth(a.x) * csch(a.x)) * a.σ)

Base.Math.cot(a::Value) = Value(cot(a.x), csc(a.x)^2 * a.σ)
Base.Math.cotd(a::Value) = Value(cotd(a.x), cscd(a.x)^2 * a.σ)
Base.Math.coth(a::Value) = Value(coth(a.x), csch(a.x)^2 * a.σ)

Base.Math.asec(a::Value) = Value(asec(a.x), a.σ / sqrt(a.x^2 - 1))
Base.Math.asecd(a::Value) = Value(asecd(a.x), a.σ / sqrt(a.x^2 - 1))
Base.Math.asech(a::Value) = Value(asech(a.x), a.σ / sqrt(1 - a.x^2))

Base.Math.acsc(a::Value) = Value(acsc(a.x), a.σ / sqrt(a.x^2 - 1))
Base.Math.acscd(a::Value) = Value(acscd(a.x), a.σ / sqrt(a.x^2 - 1))
Base.Math.acsch(a::Value) = Value(acsch(a.x), a.σ / sqrt(a.x^2 + 1))

Base.Math.acot(a::Value) = Value(acot(a.x), a.σ / (a.x^2 + 1))
Base.Math.acotd(a::Value) = Value(acotd(a.x), a.σ / (a.x^2 + 1))
Base.Math.acoth(a::Value) = Value(acoth(a.x), a.σ / abs(1 - a.x^2))

Base.Math.sinc(a::Value) = Value(sinc(a.x), a.σ * abs(cosc(a.x)))
Base.Math.cosc(a::Value) = Value(cosc(a.x), a.σ * abs((2 * a.x * cos(a.x) + (a.x^2 - 2) * sin(a.x)) / a.x^3))

Base.Math.deg2rad(a::Value) = Value(deg2rad(a.x), deg2rad(a.σ))
Base.Math.rad2deg(a::Value) = Value(rad2deg(a.x), rad2deg(a.σ))

Base.exp(a::Value) = Value(exp(a.x), exp(a.x) * a.σ)
Base.exp2(a::Value) = Value(exp2(a.x), exp2(a.x) * a.σ * log(2))
Base.exp10(a::Value) = Value(exp10(a.x), exp10(a.x) * a.σ * log(10))

Base.log(a::Value) = Value(log(a.x), a.σ / abs(a.x))
Base.log(b::T, a::Value) where {T <: Real} = Value(log(b, a.x), a.σ / abs(a.x * log(b)))
Base.log2(a::Value) = Value(log2(a.x), a.σ / abs(a.x * log(2)))
Base.log10(a::Value) = Value(log10(a.x), a.σ / abs(a.x * log(10)))
Base.log1p(a::Value) = Value(log1p(a.x), a.σ / abs(1 + a.x))

Base.sqrt(a::Value) = Value(sqrt(a.x), a.σ / (2 * sqrt(a.x)))
Base.hypot(a::Value, b::T) where {T <: Real} = Value(hypot(a.x, b), abs(a.x) * a.σ / hypot(a.x, b))
Base.hypot(a::T, b::Value) where {T <: Real} = Value(hypot(a, b.x), abs(b.x) * b.σ / hypot(a, b.x))
function Base.hypot(a::Value, b::Value)
    h = hypot(a.x, b.x)
    Value(h, hypot(a.x * a.σ, b.x * b.σ) / h)
end

Base.cbrt(a::Value) = Value(cbrt(a.x), a.σ / (3 * cbrt(a.x)^2))

## OTHER FUNCTIONS

Base.sign(a::Value) = sign(a.x)
Base.signbit(a::Value) = signbit(a.x)

Base.inv(a::Value) = Value(inv(a.x), a.σ * inv(a.x)^2)

Base.fld(a::Value, b::T) where {T <: Real} = Value(fld(a.x, b), fld(a.σ, b))
Base.fld(a::T, b::Value) where {T <: Real} = Value(fld(a, b.x), fld(b.σ * a, b.x^2))
Base.fld(a::Value, b::Value) = Value(fld(a.x, b.x), floor(typeof(a.x), hypot(fld(a.σ, b.x), fld(b.σ * a.x, b.x^2))))

Base.cld(a::Value, b::T) where {T <: Real} = Value(cld(a.x, b), cld(a.σ, b))
Base.cld(a::T, b::Value) where {T <: Real} = Value(cld(a, b.x), cld(b.σ * a, b.x^2))
Base.cld(a::Value, b::Value) = Value(cld(a.x, b.x), ceil(typeof(a.x), hypot(cld(a.σ, b.x), cld(b.σ * a.x, b.x^2))))

Base.cmp(a::Value, b::T) where {T <: Real} = cmp(a.x, b)
Base.cmp(a::T, b::Value) where {T <: Real} = cmp(a, b.x)
Base.cmp(a::Value, b::Value) = cmp(a.x, b.x)

Base.isless(a::Value, b::T) where {T <: Real} = cmp(a.x, b) < 0
Base.isless(a::T, b::Value) where {T <: Real} = cmp(a, b.x) < 0
Base.isless(a::Value, b::Value) = cmp(a.x, b.x) < 0

function Base.clamp(x::T, lh::Value) where {T <: Real}
    if x > lh.x + lh.σ
        return lh.x + lh.σ
    elseif x < lh.x - lh.σ
        return lh.x - lh.σ
    end

    x
end

function Base.clamp(rh::Value, lh::Value) where {T <: Real}
    if rh.x > lh.x + lh.σ
        return Value(lh.x + lh.σ, rh.σ)
    elseif rh.x < lh.x - lh.σ
        return Value(lh.x - lh.σ, rh.σ)
    end

    rh
end

#check how to return also the uncertainty!
Base.max(a::Value) = a.x + a.σ
Base.max(a::Value, b::T) where {T <: Real} = cmp(a.x, b) <= 0 ? b : a
Base.max(a::T, b::Value) where {T <: Real} = cmp(a, b.x) <= 0 ? b : a
Base.max(a::Value, b::Value)= cmp(a.x, b.x) <= 0 ? b : a

Base.min(a::Value) = a.x - a.σ
Base.min(a::Value, b::T) where {T <: Real} = cmp(a.x, b) <= 0 ? a : b
Base.min(a::T, b::Value) where {T <: Real} = cmp(a, b.x) <= 0 ? a : b
Base.min(a::Value, b::Value)= cmp(a.x, b.x) <= 0 ? a : b

Base.isapprox(a::Value, b::T; significance::Int = 1) where {T <: Real} =
    a.x - significance * a.σ <= b <= a.x + significance * a.σ
Base.isapprox(a::T, b::Value; significance::Int = 1) where {T <: Real} =
    b.x - significance * b.σ <= a <= b.x + significance * b.σ
Base.isapprox(a::Value, b::Value; significance::Int = 1) =
    (b.x - significance * b.σ <= a.x <= b.x + significance * b.σ) |
    (a.x - significance * a.σ <= b.x <= a.x + significance * a.σ)

Base.length(a::Value) = length(a.x)

Base.iterate(a::Value) = (a, nothing)
Base.iterate(a::Value, b::Nothing) = (a, 2)

## PRINTED RESULTS# https://docs.julialang.org/en/v1/base/math/
# precision(a::T) where {T <: Real} = floor(Int, log10(a)) - 1
Base.precision(a::Value) = floor(Int, log10(a.σ)) - 1

function Base.show(io::IO, a::Value)
    p = precision(a)
    fix_x = round(a.x * 10.0^(-p); digits = 2, base = 10) / 10
    fix_ux = round(a.σ * 10.0^(-p); digits = 2, base = 10) / 10
    p += 1
    if p == 0
        @printf(io, "%.1f ± %.1f\n", fix_x, fix_ux)
    else
        @printf(io, "(%.1f ± %.1f)·10^%d\n", fix_x, fix_ux, p)
    end
end

## NEW FUNCTIONS
val(item::Number) = item
val(item::Value) = item.x
function val(A::Array)
    type = typeof(A[1])
    if type <: Value
        type = typeof(A[1].x)
    end

    B = Array{type, ndims(A)}(undef, size(A))
    for index in eachindex(A)
        B[index] = val(A[index])
    end

    B
end

unc(item::Number) = zero(item)
unc(item::Value) = item.σ
function unc(A::Array)
    type = typeof(A[1])
    if type <: Value
        type = typeof(A[1].x)
    end

    B = Array{type, ndims(A)}(undef, size(A))
    for index in eachindex(A)
        B[index] = unc(A[index])
    end

    B
end

function set_unc(A::Array{T, N}, unc::Real) where {T, N}
    type = typeof(A[1])
    if type <: Value
        type = typeof(A[1].x)
    end

    B = Array{Value{type}, ndims(A)}(undef, size(A))
    for index in eachindex(A)
        B[index] = Value(val(A[index]), unc)
    end

    B
end

function set_unc(A::Array{T, N}, unc::Array{T2, N}) where {T, T2 <: Real, N}
    type = typeof(A[1])
    if type <: Value
        type = typeof(A[1].x)
    end

    B = Array{Value{type}, ndims(A)}(undef, size(A))
    for index in eachindex(A)
        B[index] = Value(val(A[index]), unc[index])
    end

    B
end

end
