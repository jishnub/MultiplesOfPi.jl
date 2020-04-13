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
	function PiExpTimes{N,T}(x) where {N,T<:Real}
		@assert N isa Integer
		new{N,T}(T(x))
	end
end

# Type T is inferred from x
PiExpTimes{N}(x::T) where {N,T<:Real} = PiExpTimes{N,T}(x)
PiExpTimes{0}(x::Real) = x
PiExpTimes{0}(p::PiExpTimes) = p

PiExpTimes{N}(p::PiExpTimes{M}) where {N,M} = PiExpTimes{N+M}(p.x)
PiExpTimes{N,T}(p::PiExpTimes{N,T}) where {N,T<:Real} = PiExpTimes{2N}(p.x)
PiExpTimes{N,T}(p::PiExpTimes{N}) where {N,T<:Real} = PiExpTimes{2N}(T(p.x))

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
Base.:(==)(p::PiExpTimes,y::Real) = float(p) == y
Base.:(==)(y::Real,p::PiExpTimes) = float(p) == y

Base.:(==)(::Irrational{:π},p::PiExpTimes) = false
Base.:(==)(p::PiExpTimes,::Irrational{:π}) = false

Base.:(==)(p::PiTimes,::Irrational{:π}) = isone(p.x)
Base.:(==)(::Irrational{:π},p::PiTimes) = isone(p.x)

for T in (:BigFloat,:Float64,:Float32,:Float16)
	@eval begin 
		function Base.$T(p::PiExpTimes{N}) where {N}
			if N >= 0
				y = foldl(*,(pi for i=1:N),init=one($T))
				return y*$T(p.x)
			else
				y = foldl(/,(pi for i=1:-N),init=one($T))
				return y*$T(p.x)
			end 
		end
		Base.$T(p::PiExpTimes{0}) = $T(p.x)
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
Base.zero(::Type{PiExpTimes{N,T}}) where {N,T} = PiExpTimes{N}(zero(T))

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

# The choices for conversion aren't unique as π^n can not be computed
# for negative integers. Here we resort to 
function Base.convert(::Type{PiExpTimes{N}},x::Real) where {N}
	den = N < 0 ? π^float(N) : π^N
	PiExpTimes{N}(x/den)
end
function Base.convert(::Type{PiExpTimes{N,T}},x::Real) where {N,T<:Real}
	den = N < 0 ? π^float(N) : π^N
	PiExpTimes{N}(T(x/den))
end
function Base.convert(::Type{PiExpTimes{N}},::Irrational{:π}) where {N}
	den = N < 1 ? π^float(N-1) : π^(N-1)
	PiExpTimes{N}(1/den)
end
function Base.convert(::Type{PiExpTimes{N,T}},::Irrational{:π}) where {N,T<:Real}
	den = N < 1 ? π^float(N-1) : π^(N-1)
	PiExpTimes{N}(T(1/den))
end
function Base.convert(::Type{PiExpTimes{N}},p::PiExpTimes{M}) where {M,N}
	den = N < M ? π^float(N-M) : π^(N-M)
	PiExpTimes{N}(p.x/den)
end
function Base.convert(::Type{PiExpTimes{N,T}},p::PiExpTimes{M}) where {M,N,T<:Real}
	den = N < M ? π^float(N-M) : π^(N-M)
	PiExpTimes{N}(T(p.x/den))
end
Base.convert(::Type{PiExpTimes{N}},p::PiExpTimes{N}) where {N} = p
function Base.convert(::Type{PiExpTimes{N,T}},p::PiExpTimes{N}) where {N,T<:Real}
	PiExpTimes{N}(T(p.x))
end
Base.convert(::Type{PiExpTimes{N,T}},p::PiExpTimes{N,T}) where {N,T<:Real} = p

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
