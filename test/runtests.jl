using Test
using Statistics

using UncValue

a = Value(3.1415, 0.0012)
b = Value(2.7182818, 3.4e-6)
c = Value(36458.246, 25.64)

@testset "constructor" begin
    @test a.x == 3.1415
    @test a.σ == 0.0012

    @test_throws ErrorException Value(4.1, -2.4)

    d = Value(5.248)
    @test d.x == 5.248
    @test d.σ == 0
end

@testset "precision" begin
    @test precision(a) == -4
    @test precision(b) == -7
    @test precision(c) == 0
end

@testset "addition" begin
    @test typeof(a + 2) <: Value
    @test typeof(2 + a) <: Value
    s = a + b
    @test typeof(s) <: Value
    @test s.x == (a.x + b.x)
    @test s.σ ≈ hypot(a.σ, b.σ)

    n = -a
    @test typeof(n) <: Value
    @test n.x == (-a.x)
    @test n.σ == a.σ
end

@testset "substraction" begin
    @test typeof(a - 2) <: Value
    @test typeof(2 - a) <: Value
    s = a - b
    @test typeof(s) <: Value
    @test s.x ≈ (a.x - b.x)
    @test s.σ ≈ hypot(a.σ, b.σ)
end

@testset "product" begin
    @test typeof(a * 2) <: Value
    @test typeof(2 * a) <: Value
    s = a * b
    @test typeof(s) <: Value
    @test s.x ≈ (a.x * b.x)
    @test s.σ ≈ hypot(a.σ * b.x, a.x * b.σ)
end

@testset "division" begin
    @test typeof(a / 2) <: Value
    @test typeof(2 / a) <: Value
    s = a / b
    @test typeof(s) <: Value
    @test s.x ≈ (a.x / b.x)
    @test s.σ ≈ hypot(a.σ / b.x, a.x * b.σ / b.x^2)
end

@testset "left division" begin
    s = a \ 2
    @test typeof(s) <: Value
    @test s.x ≈ (a.x \ 2)
    @test s.σ ≈ a.x^2 \ 2 * a.σ

    s = 2 \ a
    @test typeof(s) <: Value
    @test s.x ≈ (2 \ a.x)
    @test s.σ ≈ 2 \ a.σ

    s = a \ b
    @test typeof(s) <: Value
    @test s.x ≈ (a.x \ b.x)
    @test s.σ ≈ hypot(a.x \ b.σ, a.x^2 \ b.x * a.σ)
end

@testset "integer division" begin
    j = Value{Int32}(2, 1)
    k = Value{Int32}(7, 2)

    s = j ÷ k.x
    @test typeof(s) <: Value
    @test s.x ≈ j.x ÷ k.x
    @test s.σ ≈ j.σ ÷ k.x

    s = k.x ÷ j
    @test typeof(s) <: Value
    @test s.x ≈ k.x ÷ j.x
    @test s.σ ≈ k.x * j.σ ÷ j.x^2

    s = k ÷ j
    @test typeof(s) <: Value
    @test s.x ≈ k.x ÷ j.x
    @test s.σ ≈ hypot(k.x * j.σ ÷ j.x^2, j.σ ÷ k.x)

    s = fld(j, k.x)
    @test typeof(s) <: Value
    @test s.x ≈ fld(j.x, k.x)
    @test s.σ ≈ fld(j.σ, k.x)

    s = fld(k.x, j)
    @test typeof(s) <: Value
    @test s.x ≈ fld(k.x, j.x)
    @test s.σ ≈ fld(k.x * j.σ, j.x^2)

    s = fld(k, j)
    @test typeof(s) <: Value
    @test s.x ≈ fld(k.x, j.x)
    @test s.σ ≈ floor(hypot(fld(k.x * j.σ, j.x^2), fld(j.σ, k.x)))

    s = cld(j, k.x)
    @test typeof(s) <: Value
    @test s.x ≈ cld(j.x, k.x)
    @test s.σ ≈ cld(j.σ, k.x)

    s = cld(k.x, j)
    @test typeof(s) <: Value
    @test s.x ≈ cld(k.x, j.x)
    @test s.σ ≈ cld(k.x * j.σ, j.x^2)

    s = cld(k, j)
    @test typeof(s) <: Value
    @test s.x ≈ cld(k.x, j.x)
    @test s.σ ≈ ceil(hypot(cld(k.x * j.σ, j.x^2), cld(j.σ, k.x)))
end

@testset "power" begin
    @test typeof(a^3) <: Value
    @test typeof(3^a) <: Value
    s = a^b
    @test typeof(s) <: Value
    @test s.x ≈ a.x^b.x
    @test s.σ ≈ hypot(b.x * a.x^(b.x-1) * a.σ, a.x^b.x * log(abs(a.x)) * b.σ)
end

@testset "equality" begin
    @test a == 3.1415
    @test (a == 3.141) == false
    @test a != 3.141

    @test 3.1415 == a
    @test (3.141 == a) == false
    @test 3.141 != a

    @test (a == b) == false
    @test a != b
end

@testset "inequality" begin
    @test a < 3.1416
    @test (a < 3.141) == false

    @test 3.1416 > a
    @test (3.141 > a) == false

    @test a <= 3.1415
    @test a <= 3.1416
    @test (a <= 3.141) == false

    @test 3.1416 >= a
    @test 3.1415 >= a
    @test (3.141 >= a) == false

    @test a > 3.141
    @test (a > 3.1416) == false

    @test 3.141 < a
    @test (3.1416 < a) == false

    @test a >= 3.1415
    @test a >= 3.1414
    @test (a >= 3.1416) == false

    @test 3.1414 <= a
    @test 3.1415 <= a
    @test (3.146 <= a) == false

    @test a > b
    @test b < a

    @test a >= a
    @test a >= b
    @test (b >= a) == false
    @test a <= a
    @test b <= a
    @test (a <= b) == false
end

@testset "abs" begin
    s = abs(a)
    @test typeof(s) <: Value
    @test s.x == a.x
    @test s.σ == a.σ

    s = abs(-a)
    @test s.x == a.x
    @test s.σ == a.σ
end

@testset "abs2" begin
    s = abs2(a)
    @test typeof(s) <: Value
    @test s.x == a.x^2
    @test s.σ == 2 * a.x * a.σ

    s = abs2(-a)
    @test s.x == a.x^2
    @test s.σ == 2 * a.x * a.σ
end

@testset "sin/cos/tan" begin
    s = sin(b)
    @test typeof(s) <: Value
    @test s.x ≈ sin(b.x)
    @test s.σ ≈ abs(cos(b.x)) * b.σ

    c = cos(b)
    @test typeof(c) <: Value
    @test c.x ≈ cos(b.x)
    @test c.σ ≈ abs(sin(b.x)) * b.σ

    t = tan(b)
    @test typeof(t) <: Value
    @test t.x ≈ tan(b.x)
    @test t.σ ≈ sec(b.x)^2 * b.σ
end

@testset "sind/cosd/tand" begin
    s = sind(b)
    @test typeof(s) <: Value
    @test s.x ≈ sind(b.x)
    @test s.σ ≈ abs(cosd(b.x)) * b.σ

    c = cosd(b)
    @test typeof(c) <: Value
    @test c.x ≈ cosd(b.x)
    @test c.σ ≈ abs(sind(b.x)) * b.σ

    t = tand(b)
    @test typeof(t) <: Value
    @test t.x ≈ tand(b.x)
    @test t.σ ≈ secd(b.x)^2 * b.σ
end

@testset "sinh/cosh/tanh" begin
    s = sinh(b)
    @test typeof(s) <: Value
    @test s.x ≈ sinh(b.x)
    @test s.σ ≈ cosh(b.x) * b.σ

    c = cosh(b)
    @test typeof(c) <: Value
    @test c.x ≈ cosh(b.x)
    @test c.σ ≈ abs(sinh(b.x)) * b.σ

    t = tanh(b)
    @test typeof(t) <: Value
    @test t.x ≈ tanh(b.x)
    @test t.σ ≈ sech(b.x)^2 * b.σ
end

@testset "sinpi/cospi" begin
    s = sinpi(b)
    @test typeof(s) <: Value
    @test s.x ≈ sinpi(b.x)
    @test s.σ ≈ π * abs(cospi(b.x)) * b.σ

    c = cospi(b)
    @test typeof(c) <: Value
    @test c.x ≈ cospi(b.x)
    @test c.σ ≈ π * abs(sinpi(b.x)) * b.σ
end

@testset "asin/acos/atan" begin
    d = Value(0.2435, 0.0658)

    s = asin(d)
    @test typeof(s) <: Value
    @test s.x ≈ asin(d.x)
    @test s.σ ≈ d.σ / sqrt(1 - d.x^2)

    c = acos(d)
    @test typeof(c) <: Value
    @test c.x ≈ acos(d.x)
    @test c.σ ≈ d.σ / sqrt(1 - d.x^2)

    t = atan(d)
    @test typeof(t) <: Value
    @test t.x ≈ atan(d.x)
    @test t.σ ≈ d.σ / (1 + d.x^2)

    t = atan(a, b)
    @test t.x ≈ atan(a.x, b.x)
    @test t.σ ≈ hypot(b.x * a.σ, a.x * b.σ) / hypot(a.x, b.x)^2

    t = atan(a, 0.1)
    @test t.x ≈ atan(a.x, 0.1)
    @test t.σ ≈ 0.1 * a.σ / hypot(a.x, 0.1)^2

    t = atan(0.1, a)
    @test t.x ≈ atan(0.1, a.x)
    @test t.σ ≈ 0.1 * a.σ / hypot(a.x, 0.1)^2
end

@testset "asind/acosd/atand" begin
    d = Value(0.2435, 0.0658)

    s = asind(d)
    @test typeof(s) <: Value
    @test s.x ≈ asind(d.x)
    @test s.σ ≈ d.σ / sqrt(1 - d.x^2)

    c = acosd(d)
    @test typeof(c) <: Value
    @test c.x ≈ acosd(d.x)
    @test c.σ ≈ d.σ / sqrt(1 - d.x^2)

    t = atand(d)
    @test typeof(t) <: Value
    @test t.x ≈ atand(d.x)
    @test t.σ ≈ d.σ / (1 + d.x^2)

    t = atand(a, b)
    @test t.x ≈ atand(a.x, b.x)
    @test t.σ ≈ hypot(b.x * a.σ, a.x * b.σ) / hypot(a.x, b.x)^2

    t = atand(a, 0.1)
    @test t.x ≈ atand(a.x, 0.1)
    @test t.σ ≈ 0.1 * a.σ / hypot(a.x, 0.1)^2

    t = atand(0.1, a)
    @test t.x ≈ atand(0.1, a.x)
    @test t.σ ≈ 0.1 * a.σ / hypot(a.x, 0.1)^2
end

@testset "asinh/acosh/atanh" begin
    d = Value(1.2435, 0.0658)

    s = asinh(d)
    @test typeof(s) <: Value
    @test s.x ≈ asinh(d.x)
    @test s.σ ≈ d.σ / sqrt(1 + d.x^2)

    c = acosh(d)
    @test typeof(c) <: Value
    @test c.x ≈ acosh(d.x)
    @test c.σ ≈ d.σ / sqrt(d.x^2 - 1)

    d -= 1
    t = atanh(d)
    @test typeof(t) <: Value
    @test t.x ≈ atanh(d.x)
    @test t.σ ≈ d.σ / (1 - d.x^2)
end

@testset "csc/sec/cot" begin
    s = csc(b)
    @test typeof(s) <: Value
    @test s.x ≈ csc(b.x)
    @test s.σ ≈ abs(cot(b.x) * csc(b.x)) * b.σ

    c = sec(b)
    @test typeof(c) <: Value
    @test c.x ≈ sec(b.x)
    @test c.σ ≈ abs(tan(b.x) * sec(b.x)) * b.σ

    t = cot(b)
    @test typeof(t) <: Value
    @test t.x ≈ cot(b.x)
    @test t.σ ≈ csc(b.x)^2 * b.σ
end

@testset "cscd/secd/cotd" begin
    s = cscd(b)
    @test typeof(s) <: Value
    @test s.x ≈ cscd(b.x)
    @test s.σ ≈ abs(cotd(b.x) * cscd(b.x)) * b.σ

    c = secd(b)
    @test typeof(c) <: Value
    @test c.x ≈ secd(b.x)
    @test c.σ ≈ abs(tand(b.x) * secd(b.x)) * b.σ

    t = cotd(b)
    @test typeof(t) <: Value
    @test t.x ≈ cotd(b.x)
    @test t.σ ≈ cscd(b.x)^2 * b.σ
end

@testset "csch/sech/coth" begin
    s = csch(b)
    @test typeof(s) <: Value
    @test s.x ≈ csch(b.x)
    @test s.σ ≈ abs(coth(b.x) * csch(b.x)) * b.σ

    c = sech(b)
    @test typeof(c) <: Value
    @test c.x ≈ sech(b.x)
    @test c.σ ≈ abs(tanh(b.x) * sech(b.x)) * b.σ

    t = coth(b)
    @test typeof(t) <: Value
    @test t.x ≈ coth(b.x)
    @test t.σ ≈ csch(b.x)^2 * b.σ
end

@testset "acsc/asec/acot" begin
    s = acsc(b)
    @test typeof(s) <: Value
    @test s.x ≈ acsc(b.x)
    @test s.σ ≈ b.σ / sqrt(b.x^2 - 1)

    c = asec(b)
    @test typeof(c) <: Value
    @test c.x ≈ asec(b.x)
    @test c.σ ≈ b.σ / sqrt(b.x^2 - 1)

    t = acot(b)
    @test typeof(t) <: Value
    @test t.x ≈ acot(b.x)
    @test t.σ ≈ b.σ / (1 + b.x^2)
end

@testset "acscd/asecd/acotd" begin
    s = acscd(b)
    @test typeof(s) <: Value
    @test s.x ≈ acscd(b.x)
    @test s.σ ≈ b.σ / sqrt(b.x^2 - 1)

    c = asecd(b)
    @test typeof(c) <: Value
    @test c.x ≈ asecd(b.x)
    @test c.σ ≈ b.σ / sqrt(b.x^2 - 1)

    t = acotd(b)
    @test typeof(t) <: Value
    @test t.x ≈ acotd(b.x)
    @test t.σ ≈ b.σ / (1 + b.x^2)
end

@testset "acsch/asech/acoth" begin
    d = Value(0.12435, 0.0658)
    s = acsch(b)
    @test typeof(s) <: Value
    @test s.x ≈ acsch(b.x)
    @test s.σ ≈ b.σ / sqrt(b.x^2 + 1)

    c = asech(d)
    @test typeof(c) <: Value
    @test c.x ≈ asech(d.x)
    @test c.σ ≈ d.σ / sqrt(1 - d.x^2)

    t = acoth(b)
    @test typeof(t) <: Value
    @test t.x ≈ acoth(b.x)
    @test t.σ ≈ b.σ / abs(1 - b.x^2)
end

@testset "deg2rad" begin
    d = deg2rad(a)
    @test typeof(d) <: Value
    @test d.x == deg2rad(a.x)
    @test d.σ == deg2rad(a.σ)

    d = rad2deg(a)
    @test typeof(d) <: Value
    @test d.x == rad2deg(a.x)
    @test d.σ == rad2deg(a.σ)
end

@testset "exp" begin
    r = exp(a)
    @test typeof(r) <: Value
    @test r.x == exp(a.x)
    @test r.σ ≈ r.x * a.σ

    r = exp2(a)
    @test typeof(r) <: Value
    @test r.x == exp2(a.x)
    @test r.σ ≈ r.x * a.σ * log(2)

    r = exp10(a)
    @test typeof(r) <: Value
    @test r.x == exp10(a.x)
    @test r.σ ≈ r.x * a.σ * log(10)
end

@testset "log" begin
    r = log(a)
    @test typeof(r) <: Value
    @test r.x == log(a.x)
    @test r.σ ≈ a.σ / a.x

    r = log2(a)
    @test typeof(r) <: Value
    @test r.x == log2(a.x)
    @test r.σ ≈ a.σ / (a.x * log(2))

    r = log10(a)
    @test typeof(r) <: Value
    @test r.x == log10(a.x)
    @test r.σ ≈ a.σ / (a.x * log(10))

    r = log(3, a)
    @test typeof(r) <: Value
    @test r.x == log(3, a.x)
    @test r.σ ≈ a.σ / (a.x * log(3))

    r = log1p(a)
    @test typeof(r) <: Value
    @test r.x == log1p(a.x)
    @test r.σ ≈ a.σ / (a.x + 1)
end

@testset "roots" begin
    r = sqrt(a)
    @test typeof(r) <: Value
    @test r.x == sqrt(a.x)
    @test r.σ == a.σ / (2 * sqrt(a.x))

    r = cbrt(a)
    @test typeof(r) <: Value
    @test r.x == cbrt(a.x)
    @test r.σ == a.σ / (3 * cbrt(a.x)^2)

    h = hypot(a, 3)
    @test typeof(h) <: Value
    @test h.x == hypot(a.x, 3)
    @test h.σ == a.σ * a.x / h

    h = hypot(3, a)
    @test typeof(h) <: Value
    @test h.x == hypot(3, a.x)
    @test h.σ == a.σ * a.x / h

    h = hypot(a, b)
    @test typeof(h) <: Value
    @test h.x == hypot(a.x, b.x)
    @test h.σ == hypot(a.σ * a.x, b.σ * b.x) / h
end

@testset "sign" begin
    @test sign(a) == sign(a.x)
    @test sign(Value(-2.4, 2.4)) == sign(-2.4)

    @test signbit(a) == signbit(a.x)
    @test signbit(Value(-2.4, 2.4)) == signbit(-2.4)
end

@testset "inv" begin
    r = inv(a)
    @test typeof(r) <: Value
    @test r.x == inv(a.x)
    @test r.σ == a.σ / a.x^2
end

@testset "approx" begin
    @test a ≈ 3.141
    @test (a ≈ 3.145) == false
    @test isapprox(a, 3.145; significance=3)

    @test 3.141 ≈ a
    @test (3.145 ≈ a) == false
    @test isapprox(3.145, a; significance=3)

    @test (a ≈ b) == false
end

@testset "cmp" begin
    @test cmp(a, 3) == 1
    @test cmp(3, a) == -1
    @test cmp(a, a) == 0
    @test cmp(a, b) == 1
    @test cmp(b, a) == -1
end

@testset "isless" begin
    @test isless(a, b) == false
    @test isless(b, a)
    @test isless(a, 3) == false
    @test isless(3, a)
end

@testset "clamp" begin
    @test clamp(3, a) == a.x - a.σ
    @test clamp(4, a) == a.x + a.σ
    @test clamp(3.1414, a) == 3.1414
end

@testset "min/max" begin
    @test min(a, b) == min(b, a) == b
    @test max(a, b) == max(b, a) == a
    @test min(a, 3) == min(3, a) == 3
    @test max(a, 4) == max(4, a) == 4

    @test typeof(min(a, 4)) <: Value
    @test typeof(max(a, 3)) <: Value

    @test max(a) == a.x + a.σ
    @test min(a) == a.x - a.σ

    @test min(a, b, c, -2, 10) == -2
    @test min(a, b, c, 34, 13, 9) == b
    @test max(a, b, c, 2, 10) == c
    @test max(a, 4, b, 2, 1) == 4
end

@testset "val/unc" begin
    @test val(a.x) == a.x
    @test unc(a.x) == 0

    @test val(a) == a.x
    @test unc(a) == a.σ

    A = fill(Value(a.x, a.σ), (5, 2, 6))

    @test eltype(A) <: Value
    @test eltype(val(A)) <: Real
    @test mean(val(A)) ≈ a.x
    @test mean(unc(A)) ≈ a.σ
end

@testset "set_unc" begin
    A = zeros(3, 7, 2)

    uncA = set_unc(A, 0.2)
    @test eltype(uncA) <: Value
    @test mean(unc(uncA)) ≈ 0.2

    uncA = set_unc(A, fill(0.3, (3, 7, 2)))
    @test eltype(uncA) <: Value
    @test mean(unc(uncA)) ≈ 0.3

    A = fill(Value(a.x, a.σ), (5, 2, 6))
    uncA = set_unc(A, 0.2)
    @test eltype(uncA) <: Value
    @test mean(unc(uncA)) ≈ 0.2

    uncA = set_unc(A, fill(0.3, (5, 2, 6)))
    @test eltype(uncA) <: Value
    @test mean(unc(uncA)) ≈ 0.3
end
