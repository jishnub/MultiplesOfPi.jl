module MultiplesOfPi

export PiExpTimes
export PiTimes
export Pi

"""
	PiExpTimes{n}(x)

Construct a number that behaves as `x*π^n` for a real `x` and an integer `n`

# Examples
```jldoctest
julia> PiExpTimes{2}(3)
3Pi^2
```

See also: [`PiTimes`](@ref), [`Pi`](@ref)
"""
struct PiExpTimes{N,T<:Real} <: AbstractIrrational
	x :: T
	function PiExpTimes{N,T}(x::T) where {N,T<:Real}
		@assert N isa Integer
		new{N,T}(x)
	end
	function PiExpTimes{N,PiExpTimes{M,T}}(p::PiExpTimes{M,T}) where {N,M,T<:Real}
		@assert N isa Integer
		new{N,PiExpTimes{M,T}}(p)
	end
end

function PiExpTimes{N,T}(x::Real) where {N,T}
	PiExpTimes{N,T}(convert(T,x))
end

# Type T is inferred from x
PiExpTimes{N}(x::T) where {N,T<:Real} = PiExpTimes{N,T}(x)
PiExpTimes{0}(x::Real) = x
function PiExpTimes{0,T}(x::Real) where {T<:Real}
	convert(T,x)
end
function PiExpTimes{0,PiExpTimes{M,T}}(p::PiExpTimes{M}) where {M,T<:Real}
	convert(PiExpTimes{M,T},p)
end
PiExpTimes{0}(p::PiExpTimes) = p

PiExpTimes{N}(p::PiExpTimes{M}) where {N,M} = PiExpTimes{N+M}(p.x)
function PiExpTimes{N,T}(p::PiExpTimes{M}) where {N,M,T<:Real}
	PiExpTimes{N+M}(convert(T,p.x))
end
PiExpTimes{N,T}(p::PiExpTimes{N,T}) where {N,T<:Real} = PiExpTimes{2N}(p.x)
function PiExpTimes{N,T}(p::PiExpTimes{N}) where {N,T<:Real}
	PiExpTimes{2N}(convert(T,p.x))
end

function PiExpTimes{N,PiExpTimes{M,T}}(p::PiExpTimes{M,R}) where {N,M,T<:Real,R<:Real}
	p_new = convert(PiExpTimes{M,T},p)
	PiExpTimes{N,PiExpTimes{M,T}}(p_new)
end
function PiExpTimes{N,PiExpTimes{M,Real}}(p::PiExpTimes{M,R}) where {N,M,R<:Real}
	PiExpTimes{N,PiExpTimes{M,R}}(p)
end
function PiExpTimes{N,PiExpTimes{M}}(p::PiExpTimes{M,T}) where {N,M,T<:Real}
	PiExpTimes{N,PiExpTimes{M,T}}(p)
end

PiExpTimes{N}(::Irrational{:π}) where {N} = PiExpTimes{N+1}(1)
PiExpTimes{N,Irrational{:π}}(::Irrational{:π}) where {N} = PiExpTimes{N+1}(1)

"""
	PiTimes(x)

Construct a number that behaves as `x*π` for a real `x`. `PiTimes` is an alias for
`PiExpTimes{1}`

# Examples
```jldoctest
julia> PiTimes(3)
3Pi
```

See also: [`PiExpTimes`](@ref), [`Pi`](@ref)
"""
const PiTimes{T<:Real} = PiExpTimes{1,T}

"""
	Pi

The number `PiTimes(1)`, numerically equivalent to `pi`. 
Using `Pi` instead of `pi` often produces results that avoid floating-point inaccuracies.

# Examples
```jldoctest
julia> sin(Pi)
0.0

julia> sin(pi)
1.2246467991473532e-16
```
"""
const Pi = PiTimes(1)

# Helper functions to get rid of nesting
netexponent(p::T) where {T} = netexponent(T)
netexponent(::Type{PiExpTimes{N,T}}) where {N,T<:PiExpTimes} = N + netexponent(T)
netexponent(::Type{PiExpTimes{N,T}}) where {N,T<:Real} = N
netexponent(::Type{PiExpTimes{N,<:Real}}) where {N} = N
netexponent(::Type{<:Real}) = 0

rootvaltype(::T) where {T} = rootvaltype(T)
rootvaltype(::Type{PiExpTimes{N,T}}) where {N,T<:PiExpTimes} = rootvaltype(T)
rootvaltype(::Type{PiExpTimes{N,T}}) where {N,T<:Real} = T
rootvaltype(::Type{T}) where {T<:Real} = T

rootval(p::PiExpTimes) = rootval(p.x)
rootval(x::Real) = x

function setconcreterootvaltype(::Type{PiExpTimes{N,T}},U::Type{<:Real}) where {N,T<:PiExpTimes}
	PiExpTimes{N,setconcreterootvaltype(T,U)}
end
function setconcreterootvaltype(::Type{PiExpTimes{N,T}},U::Type{<:Real}) where {N,T<:Real}
	PiExpTimes{N,T}
end
function setconcreterootvaltype(::Type{PiExpTimes{N}},U::Type{<:Real}) where {N}
	PiExpTimes{N,U}
end

function rejectUnionAllroot(T::Type{<:Real},p::Real)
	rejectUnionAllroot(T,rootvaltype(p))
end
function rejectUnionAllroot(T::Type{<:PiExpTimes},U::Type{<:Real})
	if rootvaltype(T) isa UnionAll
		msg ="UnionAll types are not allowed. "
		newtype = setconcreterootvaltype(T,U)
		msg *= "Consider converting to $newtype"
		throw(ArgumentError(msg))
	end
end
function rejectUnionAllroot(T::Type{<:Real},U::Type{<:Real})
	if rootvaltype(T) isa UnionAll
		msg ="UnionAll types are not allowed. "
		throw(ArgumentError(msg))
	end
end

simplify(p::PiExpTimes) = p
function simplify(p::PiExpTimes{<:Any,<:PiExpTimes})
	T = rootvaltype(p)
	N = netexponent(p)
	x = rootval(p)
	PiExpTimes{N,T}(x)
end

# Comparison

Base.:(<)(p::PiExpTimes{N},q::PiExpTimes{N}) where {N} = p.x < q.x
function Base.:(<)(p::PiExpTimes{M},q::PiExpTimes{N}) where {M,N}
	p.x == q.x ? M < N : float(p) < float(q)
end

Base.:(<)(p::PiTimes,::Irrational{:π}) = p.x < one(p.x)
Base.:(<)(::Irrational{:π},p::PiTimes) = p.x > one(p.x)

Base.:(==)(p::PiExpTimes{N},q::PiExpTimes{N}) where {N} = p.x == q.x
function Base.:(==)(p::PiExpTimes,q::PiExpTimes)
	iszero(p.x) && iszero(q.x) ? true : float(p) == float(q)
end
Base.:(==)(p::PiExpTimes{<:Any,<:PiExpTimes},q::PiExpTimes) = simplify(p) == q
Base.:(==)(p::PiExpTimes,q::PiExpTimes{<:Any,<:PiExpTimes}) = p == simplify(q)
function Base.:(==)(p::PiExpTimes{<:Any,<:PiExpTimes},
	q::PiExpTimes{<:Any,<:PiExpTimes})
	
	simplify(p) == simplify(q)
end

Base.:(==)(p::PiExpTimes,y::Real) = float(p) == y
Base.:(==)(y::Real,p::PiExpTimes) = float(p) == y

Base.:(==)(::Irrational{:π},p::PiExpTimes) = false
Base.:(==)(p::PiExpTimes,::Irrational{:π}) = false

Base.:(==)(p::PiTimes,::Irrational{:π}) = isone(p.x)
Base.:(==)(::Irrational{:π},p::PiTimes) = isone(p.x)

for T in (:Float64,:Float32,:Float16)
	@eval begin 
		function Base.$T(p::PiExpTimes{N}) where {N}
			if N > 0
				return $T(p.x*π^N)
			else
				return $T(p.x*inv(π)^-N)
			end
		end
		Base.$T(p::PiExpTimes{0}) = $T(p.x)
		function Base.$T(p::PiTimes)
			if isone(p.x)
				return $T(π)
			else
				return $T(p.x*π)
			end
		end
	end
end

# BigFloats need to be handled differently to gain accuracy
# Convert π to BigFloat before multiplying
function Base.BigFloat(p::PiExpTimes{N}) where {N}
	if N > 0
		return BigFloat(p.x)*BigFloat(π)^N
	else
		return BigFloat(p.x)*inv(BigFloat(π))^-N
	end
end
Base.BigFloat(p::PiExpTimes{0}) = BigFloat(p.x)
function Base.BigFloat(p::PiTimes)
	if isone(p.x)
		return BigFloat(π)
	else
		return BigFloat(p.x)*π
	end
end

for f in (:iszero,:isfinite,:isnan)
	@eval Base.$f(p::PiExpTimes) = $f(p.x)
end

# Unfortunately Irrational numbers do not have a multiplicative identity of the same type,
# so we make do with something that works
Base.one(::T) where {T<:PiExpTimes} = one(T)
Base.one(::Type{<:PiExpTimes}) = true
Base.one(::Type{<:PiExpTimes{0,T}}) where {T} = one(T)

Base.zero(::T) where {T<:PiExpTimes} = zero(T)
Base.zero(::Type{PiExpTimes{N,T}}) where {N,T} = PiExpTimes{N,T}(zero(T))

Base.sign(p::PiExpTimes) = sign(p.x)
Base.signbit(p::PiExpTimes) = signbit(p.x)

# Define trigonometric functions

@inline Base.sin(p::PiTimes) = sinpi(p.x)
@inline Base.cos(p::PiTimes) = cospi(p.x)
@inline Base.sincos(p::PiTimes) = (sin(p),cos(p))

@inline Base.cis(p::PiTimes) = Complex(cos(p),sin(p))

Base.sinc(p::PiExpTimes) = sinc(float(p))
Base.sinc(p::PiExpTimes{0}) = sinc(p.x)

function Base.tan(p::PiTimes{T}) where {T}
	iszero(p.x) && return p.x
	r = abs(rem(p.x,one(p.x)))
	iszero(r) && return copysign(zero(T),p.x)
	2r == one(r) && return copysign(T(Inf),p.x)
	4r == one(r) && return copysign(one(T),p.x)
	tan(p.x*π)
end

# Hyperbolic functions

function Base.tanh(z::Complex{<:PiTimes})
	iszero(real(z)) && return im*tan(imag(z))
	tanh(float(z))
end

# Arithmetic operators

Base.:(+)(p1::PiExpTimes{N},p2::PiExpTimes{N}) where {N} = PiExpTimes{N}(p1.x + p2.x)
Base.:(+)(p1::PiExpTimes{0},p2::PiExpTimes{0}) = p1.x + p2.x

Base.:(+)(p::PiTimes,::Irrational{:π}) = PiTimes(p.x + one(p.x))
Base.:(+)(::Irrational{:π},p::PiTimes) = PiTimes(p.x + one(p.x))

Base.:(-)(p::PiExpTimes{N}) where {N} = PiExpTimes{N}(-p.x)
Base.:(-)(p::PiExpTimes{0}) = -p.x

Base.:(-)(p1::PiExpTimes{N},p2::PiExpTimes{N}) where {N} = PiExpTimes{N}(p1.x-p2.x)
Base.:(-)(p1::PiExpTimes{0},p2::PiExpTimes{0}) = p1.x-p2.x

Base.:(-)(p::PiTimes,::Irrational{:π}) = PiTimes(p.x - one(p.x))
Base.:(-)(::Irrational{:π},p::PiTimes) = PiTimes(one(p.x) - p.x)

Base.:(/)(p1::PiExpTimes{N},p2::PiExpTimes{N}) where {N} = p1.x/p2.x
Base.:(/)(p1::PiExpTimes{M},p2::PiExpTimes{N}) where {M,N} = PiExpTimes{M-N}(p1.x/p2.x)

Base.:(/)(y::Real,p::PiExpTimes) = y*inv(p)
Base.:(/)(p1::PiExpTimes{N},y::Real) where {N} = PiExpTimes{N}(p1.x/y)

Base.:(/)(::Irrational{:π},p::PiExpTimes{N}) where {N} = PiExpTimes{1-N}(inv(p.x))
Base.:(/)(::Irrational{:π},p::PiTimes) = inv(p.x)
Base.:(/)(p::PiExpTimes{N},::Irrational{:π}) where {N} = PiExpTimes{N-1}(p.x)
Base.:(/)(p::PiTimes,::Irrational{:π})= p.x

function Base.:(/)(p::PiExpTimes{N},z::Complex{<:PiExpTimes{M}}) where {M,N}
	are,aim = p.x, zero(p.x)
	bre,bim = real(z).x, imag(z).x
	Complex(are,aim)/Complex(bre,bim) * PiExpTimes{N-M}(1)
end

function Base.:(/)(p::PiExpTimes{N},z::Complex{<:PiExpTimes{N}}) where {N}
	are,aim = p.x, zero(p.x)
	bre,bim = real(z).x, imag(z).x
	Complex(are,aim)/Complex(bre,bim)
end

function Base.:(/)(a::Complex{<:PiExpTimes{N}},b::Complex{<:Real}) where {N}
	are,aim = real(a).x, imag(a).x
	bre,bim = reim(b)
	Complex(are,aim)/Complex(bre,bim) * PiExpTimes{N}(1)
end

function Base.:(/)(a::Complex{<:PiExpTimes{N}},b::Complex{<:PiExpTimes{M}}) where {M,N}
	are,aim = real(a).x, imag(a).x
	bre,bim = real(b).x, imag(b).x
	Complex(are,aim)/Complex(bre,bim) * PiExpTimes{N-M}(1)
end

function Base.:(/)(a::Complex{<:PiExpTimes{N}},b::Complex{<:PiExpTimes{N}}) where {N}
	are,aim = real(a).x, imag(a).x
	bre,bim = real(b).x, imag(b).x
	Complex(are,aim)/Complex(bre,bim)
end

Base.inv(p::PiExpTimes{N}) where {N} = PiExpTimes{-N}(inv(p.x))

Base.:(*)(p1::PiExpTimes{N},y::Real) where {N} = PiExpTimes{N}(p1.x*y)
Base.:(*)(y::Real,p1::PiExpTimes{N}) where {N} = PiExpTimes{N}(p1.x*y)

Base.:(*)(b::Bool,p::PiExpTimes) = ifelse(b,p,flipsign(zero(p),p.x))
Base.:(*)(p::PiExpTimes,b::Bool) = ifelse(b,p,flipsign(zero(p),p.x))

function Base.:(*)(p1::PiExpTimes{M},p2::PiExpTimes{N}) where {M,N}
	PiExpTimes{M+N}(p1.x*p2.x)
end

Base.:(*)(z::Complex{Bool},p::PiExpTimes) = Complex(real(z)*p,imag(z)*p)
Base.:(*)(p::PiExpTimes,z::Complex{Bool}) = Complex(real(z)*p,imag(z)*p)

Base.:(*)(::Irrational{:π},p::PiExpTimes) = PiTimes(p)
Base.:(*)(p::PiExpTimes,::Irrational{:π}) = PiTimes(p)

# Not type-stable!
Base.:(^)(p::PiExpTimes{N},n::Integer) where {N} = PiExpTimes{N*n}(p.x^n)

for op in Symbol[:/,:*]
	@eval Base.$op(p::PiExpTimes,x::AbstractIrrational) = $op(Float64(p),Float64(x))
	@eval Base.$op(x::AbstractIrrational,p::PiExpTimes) = $op(Float64(x),Float64(p))
end

# Rational divide
Base.:(//)(p::PiExpTimes{N}, n::Real) where {N} = PiExpTimes{N}(p.x//n)
Base.:(//)(n::Real, p::PiExpTimes{N}) where {N} = PiExpTimes{-N}(n//p.x)
Base.:(//)(p::PiExpTimes{M},q::PiExpTimes{N}) where {M,N} = PiExpTimes{M-N}(p.x//q.x)
Base.:(//)(p::PiExpTimes{N},q::PiExpTimes{N}) where {N} = p.x//q.x

Base.:(//)(::Irrational{:π},p::PiExpTimes{N}) where {N} = PiExpTimes{1-N}(1//p.x)
Base.:(//)(p::PiExpTimes{N},::Irrational{:π}) where {N} = PiExpTimes{N-1}(p.x//1)
Base.:(//)(::Irrational{:π},p::PiTimes) = 1//p.x
Base.:(//)(p::PiTimes,::Irrational{:π}) = p.x//1

Base.:(//)(p::Complex{<:PiExpTimes},q::PiExpTimes) = Complex(real(p)//q,imag(p)//q)

function Base.:(//)(p::PiExpTimes{N},z::Complex{<:PiExpTimes{M}}) where {M,N}
	are,aim = p.x, zero(p.x)
	bre,bim = real(z).x, imag(z).x
	Complex(are,aim)//Complex(bre,bim) * PiExpTimes{N-M}(1)
end

function Base.:(//)(p::PiExpTimes{N},z::Complex{<:PiExpTimes{N}}) where {N}
	are,aim = p.x, zero(p.x)
	bre,bim = real(z).x, imag(z).x
	Complex(are,aim)//Complex(bre,bim)
end

function Base.:(//)(a::Complex{<:PiExpTimes{N}},b::Complex{<:PiExpTimes{M}}) where {M,N}
	are,aim = real(a).x, imag(a).x
	bre,bim = real(b).x, imag(b).x
	Complex(are,aim)//Complex(bre,bim) * PiExpTimes{N-M}(1)
end

function Base.:(//)(a::Complex{<:PiExpTimes{N}},b::Complex{<:PiExpTimes{N}}) where {N}
	are,aim = real(a).x, imag(a).x
	bre,bim = real(b).x, imag(b).x
	Complex(are,aim)//Complex(bre,bim)
end

# Conversion and promotion
function Base.promote_rule(::Type{PiExpTimes{N,R}},::Type{PiExpTimes{N,S}}) where {N,R,S}
	PiExpTimes{N,promote_type(R,S)}
end

Base.promote_rule(::Type{PiTimes{T}}, ::Type{Irrational{:π}}) where {T} = PiTimes{T}
Base.promote_rule(::Type{Irrational{:π}}, ::Type{PiTimes{T}}) where {T} = PiTimes{T}

function Base.promote_rule(::Type{Complex{PiTimes{T}}}, ::Type{Irrational{:π}}) where {T}
	Complex{PiTimes{T}}
end
function Base.promote_rule(::Type{Irrational{:π}}, ::Type{Complex{PiTimes{T}}}) where {T}
	Complex{PiTimes{T}}
end

# Conversions

struct ExponentMismatchError end
Base.showerror(io::IO,e::ExponentMismatchError) = print(io,"Exponent mismatch in type")
struct IncompatibleTypesError end
Base.showerror(io::IO,e::IncompatibleTypesError) = print(io,
	"Incompatible types")

# Float to PiExpTimes or PiTimes

# PiExpTimes{0} is a special case, equivalent to identity
Base.convert(::Type{<:PiExpTimes{0}},x::Real) = x
Base.convert(::Type{PiExpTimes{0,T}},x::Real) where {T<:Real} = convert(T,x)

function Base.convert(::Type{<:PiExpTimes{N,P}},x::T) where {N,T<:Real,P<:Real}
	# eg. convert(PiExpTimes{2,Int},2) == Int(2Pi^-2)*Pi^2
	# This fails as Pi^-2 is not an Int
	throw(IncompatibleTypesError())
end
# Conversion involving Real introduces nesting
function Base.convert(::Type{PiExpTimes{N,Real}},x::T) where {N,T<:Real}
	# eg. convert(PiExpTimes{2,Real},2) == 2Pi^-2*Pi^2
	convert(PiExpTimes{N},x)
end
function Base.convert(::Type{<:PiExpTimes{N,P}},x::T) where {N,T<:Real,P<:AbstractFloat}
	# eg. convert(PiExpTimes{2,Float64},2) == Float64(2/Pi^2)*Pi^2
	PiExpTimes{N,P}(convert(P,PiExpTimes{-N}(x)))
end
function Base.convert(::Type{<:PiTimes{P}},x::T) where {T<:Real,P<:Real}
	# eg. convert(PiTimes{Int},2) == Int(2Pi^-1)*Pi
	# This fails as Pi^-1 is not an Int
	# This method is provided for ambiguity resolution
	throw(IncompatibleTypesError())
end
function Base.convert(::Type{<:PiTimes{P}},x::T) where {T<:Real,P<:AbstractFloat}
	# eg. convert(PiTimes{Float64},2) == Float64(2Pi^-1)*Pi
	PiTimes{P}(convert(P,PiExpTimes{-1}(x)))
end
function Base.convert(::Type{<:PiExpTimes{N}},x::T) where {N,T<:Real}
	# eg. convert(PiExpTimes{2},2)
	p_new = PiExpTimes{-N,T}(x)
	PiExpTimes{N,PiExpTimes{-N,T}}(p_new)
end
function Base.convert(::Type{<:PiTimes},x::T) where {T<:Real}
	# eg. convert(PiTimes,2) == 2Pi^-1*Pi
	p_new = PiExpTimes{-1,T}(x)
	PiTimes{PiExpTimes{-1,T}}(p_new)
end

# Conversion from real to nested
function Base.convert(::Type{PiExpTimes{N,PiExpTimes{M,R}}},x::T) where {N,M,T<:Real,R<:Real}
	# eg. convert(PiExpTimes{1,PiExpTimes{1,Float64}},2)
	# Error checks
	rejectUnionAllroot(PiExpTimes{N,PiExpTimes{M,R}},x)

	PiExpTimes{N,PiExpTimes{M,R}}(PiExpTimes{M,R}(x))
end
function Base.convert(::Type{PiTimes{PiExpTimes{M,R}}},x::T) where {M,T<:Real,R<:Real}
	# Error checks
	rejectUnionAllroot(PiTimes{PiExpTimes{M,R}},x)
	PiTimes{PiExpTimes{M,R}}(PiExpTimes{M,R}(x))
end

# Irrational{π} to PiExpTimes or PiTimes

# General case that infers the type
function Base.convert(::Type{<:PiExpTimes{N}},::Irrational{:π}) where {N}
	convert(PiExpTimes{N},Pi)
end

# Special case, convert to bits types
function Base.convert(::Type{PiExpTimes{N,T}},::Irrational{:π}) where {N,T<:Real}
	# eg. convert(PiExpTimes{2,Int},π) == Int(1/π)*Pi^2
	# This fails as 1/π is not an Int
	throw(IncompatibleTypesError())
end
function Base.convert(::Type{PiExpTimes{N,T}},::Irrational{:π}) where {N,T<:AbstractFloat}
	# eg. convert(PiExpTimes{2,Float64},π) == Float64(1/π)*Pi^2
	# Error checks
	rejectUnionAllroot(PiExpTimes{N,T},Pi)
	PiExpTimes{N,T}(convert(T,PiExpTimes{1-N}(1)))
end
# Conversion to Real introduces nesting
function Base.convert(::Type{PiExpTimes{N,Real}},::Irrational{:π}) where {N}
	# eg. convert(PiExpTimes{2,Real},π) == Pi^-1*Pi^2
	PiExpTimes{N,PiExpTimes{1-N,Int}}(PiExpTimes{1-N,Int}(1))
end

# Conversion to nested types from π
function Base.convert(::Type{PiExpTimes{N,PiExpTimes{M,R}}},
	::Irrational{:π}) where {N,M,R<:Real}
	# eg. convert(PiExpTimes{2,PiExpTimes{-1,Int}},π) == Pi^-1*Pi^2
	# Error checks
	rejectUnionAllroot(PiExpTimes{N,PiExpTimes{M,R}},Pi)

	p_new = convert(PiExpTimes{M,R},PiExpTimes{1-N}(1))
	PiExpTimes{N,PiExpTimes{M,R}}(p_new)
end
function Base.convert(::Type{PiExpTimes{N,PiExpTimes{M}}},
	::Irrational{:π}) where {N,M,R<:Real}

	p_new = convert(PiExpTimes{M},PiExpTimes{1-N}(1))
	convert(PiExpTimes{N,PiExpTimes{M}},p_new)
end

# Conversion to PiTimes from π
Base.convert(::Type{<:PiTimes},::Irrational{:π}) = Pi
function Base.convert(::Type{PiTimes{T}},::Irrational{:π}) where {T<:Real}
	rejectUnionAllroot(PiTimes{T},Pi)
	PiTimes{T}(convert(T,1))
end
function Base.convert(::Type{PiTimes{T}},::Irrational{:π}) where {T<:AbstractFloat}
	rejectUnionAllroot(PiTimes{T},Pi)
	PiTimes{T}(convert(T,1))
end

# Conversion from PiExpTimes to PiExpTimes
# General case without types specified, will be inferred
# Safe to run this as there won't be an error
function Base.convert(::Type{<:PiExpTimes{N}},p::PiExpTimes{M,T}) where {M,N,T<:Real}
	# eg. convert(PiExpTimes{2},Pi)
	p_new = PiExpTimes{M-N,T}(p.x)
	PiExpTimes{N,PiExpTimes{M-N,T}}(p_new)
end
function Base.convert(::Type{<:PiTimes},p::PiExpTimes{M,T}) where {M,T<:Real}
	# eg. convert(PiExpTimes{2},Pi)
	p_new = PiExpTimes{M-1,T}(p.x)
	PiTimes{PiExpTimes{M-1,T}}(p_new)
end

# General case with type and exponent conversion
function Base.convert(::Type{PiExpTimes{N,T}},p::PiExpTimes{M,R}) where {M,N,T<:Real,R<:Real}
	# eg. convert(PiExpTimes{2,Int},Pi)
	# This should fail as the result is Pi^-1*Pi^2 and Pi^-1 is not an Int
	throw(IncompatibleTypesError())
end
function Base.convert(::Type{PiExpTimes{N,Real}},p::PiExpTimes{M,R}) where {M,N,R<:Real}
	# eg. convert(PiExpTimes{2,Real},Pi) == Pi^-1*Pi^2
	convert(PiExpTimes{N},p)
end
# Special case where the type is float
# In this case the root value will be converted to float
function Base.convert(::Type{PiExpTimes{N,T}},
	p::PiExpTimes{M,R}) where {M,N,T<:AbstractFloat,R<:Real}
	# eg. convert(PiExpTimes{2,Float64},Pi) == Float64(1/Pi)*Pi^2
	# Error checks
	rejectUnionAllroot(PiExpTimes{N,T},p)
	
	p_new = PiExpTimes{M-N}(p.x)
	PiExpTimes{N}(convert(T,p_new))
end

# Special cases with same exponent
# Most general special case without type specified
Base.convert(::Type{<:PiExpTimes{N}},p::PiExpTimes{N}) where {N} = p
# Specified types must be compatible otherwise
function Base.convert(::Type{PiExpTimes{N,T}},p::PiExpTimes{N,R}) where {N,T<:Real,R<:Real}
	# Error checks
	rejectUnionAllroot(PiExpTimes{N,T},p)
	PiExpTimes{N}(convert(T,p.x))
end
Base.convert(::Type{PiExpTimes{N,T}},p::PiExpTimes{N,T}) where {N,T<:Real} = p
Base.convert(::Type{PiExpTimes{N,T}},p::PiExpTimes{N,T}) where {N,T<:AbstractFloat} = p

# General nested case
function Base.convert(::Type{PiExpTimes{N,PiExpTimes{M,R}}},
	p::PiExpTimes{Q,T}) where {N,M,Q,T<:Real,R<:Real}
	
	rejectUnionAllroot(PiExpTimes{N,PiExpTimes{M,R}},p)

	# How this works: convert(PiExpTimes{2,PiExpTimes{2,Int}},Pi^4) == Pi^2*Pi^2
	p_new = convert(PiExpTimes{M,R},PiExpTimes{Q-N,T}(p.x))
	PiExpTimes{N,PiExpTimes{M,R}}(p_new)
end
# Special case
function Base.convert(::Type{PiTimes{PiExpTimes{M,R}}},
	p::PiExpTimes{Q,T}) where {M,Q,T<:Real,R<:Real}
	
	netexponent(PiTimes{PiExpTimes{M,R}}) == netexponent(p) ||
		throw(ExponentMismatchError())
	rejectUnionAllroot(PiTimes{PiExpTimes{M,R}},p)

	# How this works: convert(PiTimes{PiExpTimes{2,Int}},Pi^3) == Pi^2*Pi
	PiTimes{PiExpTimes{M,R}}(convert(PiExpTimes{M,R},PiExpTimes{Q-1,T}(p.x)))
end
# Even more special case
function Base.convert(::Type{PiTimes{PiExpTimes{M,R}}},
	p::PiTimes) where {M,R<:Real}
	
	netexponent(PiTimes{PiExpTimes{M,R}}) == netexponent(p) ||
		throw(ExponentMismatchError())
	rejectUnionAllroot(PiTimes{PiExpTimes{M,R}},p)

	# How this works: convert(PiTimes{PiExpTimes{0,Int}},Pi) == Pi^0*Pi
	p_new = convert(PiExpTimes{M,R},p)
	PiTimes{PiExpTimes{M,R}}(p_new)
end
function Base.convert(::Type{PiTimes{PiExpTimes{M,R}}},
	p::PiTimes{PiExpTimes{M,R}}) where {M,R<:Real} 
	p
end

# Conversions to PiTimes
Base.convert(::Type{<:PiTimes},p::PiTimes) = p
Base.convert(::Type{PiTimes{T}},p::PiTimes{T}) where {T<:Real} = p
Base.convert(::Type{PiTimes{T}},p::PiTimes{T}) where {T<:AbstractFloat} = p
function Base.convert(::Type{PiTimes{T}},p::PiTimes{R}) where {T<:AbstractFloat,R<:Real}
	PiTimes{T}(convert(T,p.x))
end
Base.convert(::Type{PiTimes{T}},p::PiTimes) where {T<:Real} = PiTimes{T}(p.x)
function Base.convert(::Type{PiTimes{T}},p::PiExpTimes{N,R}) where {N,T<:Real,R<:Real}
	# This will catch cases like convert(PiTimes{Int},Pi^2) = Int(1/Pi)*Pi^2
	# This fails because 1/Pi is not an Int
	throw(IncompatibleTypesError())
end
function Base.convert(::Type{PiTimes{Real}},p::PiExpTimes{N,R}) where {N,R<:Real}
	# This will catch cases like convert(PiTimes{Real},Pi^2) = Pi^-1*Pi^2
	PiTimes{PiExpTimes{N-1,R}}(PiExpTimes{N-1,R}(p.x))
end
function Base.convert(::Type{PiTimes{T}},p::PiExpTimes{N,R}) where {N,T<:AbstractFloat,R<:Real}
	# This will catch cases like convert(PiTimes{Float64},Pi^2) = Float64(Pi)*Pi
	PiTimes{T}(convert(T,PiExpTimes{N-1,R}(p.x)))
end

function Base.show(io::IO,p::PiExpTimes{N}) where {N}
	x = p.x

	expstr(p) = isone(N) ? "Pi" : "Pi^"*string(N)
	
	tostr(x,p) = isone(x) ? expstr(p) : string(x)*"*"*expstr(p)
	tostr(x::Integer,p) = isone(x) ? expstr(p) : string(x)*expstr(p)
	tostr(x::AbstractFloat,p) = isone(x) ? expstr(p) : 
								isinf(x) || isnan(x) ? string(x) :
								string(x)*"*"*expstr(p)
	tostr(x::Rational,p) = "("*string(x)*")"*expstr(p)

	str = iszero(x) ? string(x) : tostr(x,p)
	print(io,str)
end

end # module
