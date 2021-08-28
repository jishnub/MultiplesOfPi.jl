module MultiplesOfPi

export PiExpTimes
export Pi

"""
	PiExpTimes(x, n)

Construct a number that behaves as `x*π^n` for a real `x` and an integer `n`

# Examples
```jldoctest
julia> PiExpTimes(2, 3)
2Pi^3
```

See also: [`Pi`](@ref)
"""
struct PiExpTimes{T<:Real} <: AbstractIrrational
	x :: T
	n :: Int

	PiExpTimes{T}(x::T, n::Int) where {T<:Real} = new(x, n)
	PiExpTimes{PiExpTimes{T}}(p::PiExpTimes{T}, n::Int) where {T<:Real} = new(p, n)
	PiExpTimes{PiExpTimes}(p::PiExpTimes, n::Int) = new(p, n)
	PiExpTimes{Irrational{:π}}(p::Irrational{:π}, n::Int) = new(p, n)
end

PiExpTimes{T}(x::Real, n::Int) where {T<:Real} = PiExpTimes{T}(T(x), n)
PiExpTimes{T}(::Irrational{:π}, n::Int) where {T<:Real} = PiExpTimes{T}(T(1), n + 1)
PiExpTimes{T}(p::PiExpTimes, n::Int) where {T<:Real} = PiExpTimes{T}(p.x, p.n + n)
PiExpTimes{T}(x::Real, n::Integer = 0) where {T<:Real} = PiExpTimes{T}(x, Int(n))

PiExpTimes(::Irrational{:π}, n::Integer = 0) = PiExpTimes{Int}(1, Int(n) + 1)
PiExpTimes(x::Real, n::Integer = 0) = PiExpTimes{typeof(x)}(x, Int(n))
# avoid nesting
PiExpTimes(p::PiExpTimes, n::Integer = 0) = PiExpTimes(p.x, p.n + n)

"""
	Pi

The number `PiExpTimes(1, 1)`, numerically equivalent to `pi`.
Using `Pi` instead of `pi` often produces results that avoid floating-point inaccuracies.

# Examples
```jldoctest
julia> sin(Pi)
0.0

julia> sin(π)
1.2246467991473532e-16

julia> exp(im*Pi) + 1 == 0
true

julia> exp(im*π) + 1 == 0
false

julia> (2//3)Pi + Pi == (5//3)Pi
true

julia> (2//3)π + π == (5//3)π
false
```
"""
const Pi = PiExpTimes(1, 1)

function Base.:(<)(p::PiExpTimes, q::PiExpTimes)
	p.x == q.x ? p.n < q.n :
	p.n == q.n ? p.x < q.x :
	float(p) < float(q)
end
function Base.:(==)(p::T, q::T) where {T<:PiExpTimes}
	iszero(p.x) && iszero(q.x) ? true :
	(p.n == q.n && p.x == q.x) ? true :
	false
end

# Compare nested types by simplifying and promoting to a common type
function Base.:(==)(p::PiExpTimes, q::PiExpTimes)
	p2, q2 = promote(simplify(p), simplify(q))
	p2 == q2
end

Base.:(==)(p::PiExpTimes, x::Real) = p == PiExpTimes(x, 0)
Base.:(==)(x::Real, p::PiExpTimes) = p == x

simplify(p) = p
simplify(p::PiExpTimes{<:PiExpTimes}) = simplify(p.x) * Pi^(p.n)

for f in [:(<), :(==)]
	@eval Base.$f(p::PiExpTimes, ::Irrational{:π}) = $f(p, Pi)
	@eval Base.$f(::Irrational{:π}, p::PiExpTimes) = $f(Pi, p)
end

Base.promote_rule(::Type{Irrational{:π}}, ::Type{PiExpTimes{T}}) where {T} = PiExpTimes{T}
Base.promote_rule(::Type{PiExpTimes{T}}, ::Type{Irrational{:π}}) where {T} = PiExpTimes{T}

Base.promote_rule(::Type{PiExpTimes{T}}, ::Type{Bool}) where {T} = PiExpTimes{T}
Base.promote_rule(::Type{PiExpTimes{T}}, ::Type{S}) where {T,S<:AbstractIrrational} = PiExpTimes{promote_type(S,T)}
Base.promote_rule(::Type{PiExpTimes{T}}, ::Type{S}) where {T,S<:Real} = PiExpTimes{promote_type(S,T)}
Base.promote_rule(::Type{PiExpTimes{T}}, ::Type{Float16}) where {T} = PiExpTimes{promote_type(Float16,T)}
Base.promote_rule(::Type{PiExpTimes{T}}, ::Type{Float32}) where {T} = PiExpTimes{promote_type(Float32,T)}

Base.promote_rule(::Type{PiExpTimes{T}}, ::Type{PiExpTimes{S}}) where {T,S} = PiExpTimes{promote_type(T, S)}

Base.promote_rule(::Type{PiExpTimes{T}}, ::Type{Complex{S}}) where {T<:Real,S<:Real} = Complex{PiExpTimes{promote_type(T, S)}}

for T in (:Float64,:Float32,:Float16)
	@eval begin
		function Base.$T(p::PiExpTimes)
			n = p.n
			if n == 0
				return $T(p.x)
			elseif n > 0
				return $T(p.x*π^n)
			else
				return $T(p.x*inv(float(π))^-n)
			end
		end
	end
end

# BigFloats need to be handled differently to gain accuracy
# Convert π to BigFloat before multiplying
function Base.BigFloat(p::PiExpTimes)
	n = p.n
	if n == 0
		return BigFloat(p.x)
	elseif n > 0
		return BigFloat(p.x)*BigFloat(π)^n
	else
		return BigFloat(p.x)*inv(BigFloat(π))^-n
	end
end

for f in (:iszero, :isfinite, :isnan)
	@eval Base.$f(p::PiExpTimes) = $f(p.x)
end

# Unfortunately Irrational numbers do not have a multiplicative identity of the same type,
# so we make do with something that works
Base.one(::Type{PiExpTimes}) = true
Base.isone(p::PiExpTimes) = (x = simplify(p); isone(x.x) && iszero(x.n))

Base.zero(p::PiExpTimes{T}) where {T} = PiExpTimes{T}(zero(T), p.n)
Base.zero(::Type{PiExpTimes{T}}) where {T} = PiExpTimes{T}(zero(T))

Base.sign(p::PiExpTimes) = sign(p.x)
Base.signbit(p::PiExpTimes) = signbit(p.x)

# Define trigonometric functions

# this is type-unstable, but Union splitting does its magic
@inline dividebypiexp(p::PiExpTimes, n) = p.n == n ? p.x : p.x * float(π)^(p.n-n)
@inline dividebypi(p) = dividebypiexp(p, 1)

Base.sin(p::PiExpTimes) = p.n == 0 ? sin(p.x) : sinpi(dividebypi(p))
Base.cos(p::PiExpTimes) = p.n == 0 ? cos(p.x) : cospi(dividebypi(p))
if VERSION > v"1.6.0"
	Base.sincos(p::PiExpTimes) = p.n == 0 ? sincos(p.x) : sincospi(dividebypi(p))
else
	Base.sincos(p::PiExpTimes) = (sin(p),cos(p))
end

function Base.cis(p::PiExpTimes)
	s, c = sincos(p);
	Complex(c, s)
end

@inline function Base.tan(p::PiExpTimes{T}) where {T}
	Tf = float(T)
	if simplify(p).n == 1
		r = abs(rem(p.x, one(p.x)))
		iszero(r) && return copysign(zero(Tf), p.x)
		2r == one(r) && return copysign(Tf(Inf), p.x)
		4r == one(r) && return copysign(one(Tf), p.x)
	end
	tan(float(p))
end

function Base.tanh(z::Complex{<:PiExpTimes})
	iszero(real(z)) && return im*tan(imag(z))
	tanh(float(z))
end

Base.sinc(p::PiExpTimes) = sinc(float(p))

# Arithmetic operators

function match_exponent(p1, p2)
	minexp = min(p1.n, p2.n)
	x = dividebypiexp(p1, minexp)
	y = dividebypiexp(p2, minexp)
	x, y, minexp
end

function Base.:(+)(p1::PiExpTimes, p2::PiExpTimes)
	x, y, exponent = match_exponent(p1, p2)
	PiExpTimes(x + y, exponent)
end

Base.:(-)(p::PiExpTimes) = PiExpTimes(-p.x, p.n)

function Base.:(-)(p1::PiExpTimes, p2::PiExpTimes)
	x, y, exponent = match_exponent(p1, p2)
	PiExpTimes(x - y, exponent)
end

Base.inv(p::PiExpTimes) = PiExpTimes(inv(p.x), -p.n)

function Base.:(*)(p1::PiExpTimes, p2::PiExpTimes)
	PiExpTimes(p1.x*p2.x, p1.n+p2.n)
end
function Base.:(/)(p1::PiExpTimes, p2::PiExpTimes)
	PiExpTimes(p1.x/p2.x, p1.n-p2.n)
end

for f in [:(+), :(-), :(*), :(/)]
	@eval Base.$f(p::PiExpTimes, ::Irrational{:π}) = $f(p, Pi)
	@eval Base.$f(::Irrational{:π}, p::PiExpTimes) = $f(Pi, p)
end

Base.:(*)(b::Bool, p::PiExpTimes) = ifelse(b, p, flipsign(zero(p), p.x))
Base.:(/)(b::Bool, p::PiExpTimes) = PiExpTimes(b / p.x, -p.n)

Base.:(^)(p::PiExpTimes, n::Integer) = PiExpTimes(p.x^n, n*p.n)

# Rational divide
Base.:(//)(p::PiExpTimes, q::PiExpTimes) = PiExpTimes(p.x//q.x, p.n - q.n)
Base.:(//)(p::PiExpTimes, y::Real) = PiExpTimes(p.x//y, p.n)
Base.:(//)(p::PiExpTimes, ::Irrational{:π}) = PiExpTimes(p.x//1, p.n-1)
Base.:(//)(y::Real, p::PiExpTimes) = PiExpTimes(y//p.x, -p.n)
Base.:(//)(::Irrational{:π}, p::PiExpTimes) = PiExpTimes(1//p.x, -p.n+1)

# Pretty-printing

function Base.show(io::IO, p::PiExpTimes)
	x = p.x

	expstr(p) = isone(p.n) ? "Pi" : "Pi^"*string(p.n)

	tostr(x,p) = isone(x) ? expstr(p) : string(x)*"*"*expstr(p)
	tostr(x::Integer,p) = isone(x) ? expstr(p) : string(x)*expstr(p)
	tostr(x::AbstractFloat,p) = isone(x) ? expstr(p) :
								isinf(x) || isnan(x) ? string(x) :
								string(x)*"*"*expstr(p)
	tostr(x::Rational,p) = "("*string(x)*")"*expstr(p)

	str = iszero(x) ? string(x) : tostr(x,p)
	print(io, str)
end

end # module
