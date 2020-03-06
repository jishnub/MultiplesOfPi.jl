module MultiplesOfPi

export PiTimes
export Pi

struct PiTimes{T<:Real} <: AbstractIrrational
	x :: T
end

const Pi = PiTimes(1)

Base.:(<)(p::PiTimes,q::PiTimes) = p.x < q.x
Base.:(==)(p::PiTimes,q::PiTimes) = p.x == q.x
Base.:(==)(p::PiTimes,y::Real) = float(p) == y
Base.:(==)(y::Real,p::PiTimes) = float(p) == y

Base.:(==)(p::PiTimes,::Irrational{:π}) = isone(p.x)
Base.:(==)(::Irrational{:π},p::PiTimes) = isone(p.x)

for T in (:BigFloat,:Float64,:Float32,:Float16)
	eval(quote
		Base.$T(p::PiTimes) = π*$T(p.x)
	end)
end

for f in (:iszero,:isfinite,:isnan)
	eval(quote
		Base.$f(p::PiTimes) = $f(p.x)
	end)
end

# Unfortunately Irrational numbers do not have a multiplicative identity of the same type,
# so we make do with something that works
# NOTE: This will be changed in the next minor release to one(::PiTimes) = true
Base.one(::Type{PiTimes{T}}) where {T} = one(T)

# Define trigonometric functions

@inline Base.sin(p::PiTimes) = sinpi(p.x)
@inline Base.cos(p::PiTimes) = cospi(p.x)
@inline Base.sincos(p::PiTimes) = (sin(p),cos(p))

@inline Base.cis(p::PiTimes) = Complex(cos(p),sin(p))

function Base.tan(p::PiTimes{T}) where {T}
	iszero(p.x) && return p.x
	r = abs(rem(p.x,one(p.x)))
	iszero(r) && return copysign(zero(T),p.x)
	2r == one(r) && return copysign(T(Inf),p.x)
	4r == one(r) && return copysign(one(T),p.x)
	tan(p.x*π)
end

# Hyperbolic functions

function Base.tanh(z::Complex{PiTimes{T}}) where {T}
	iszero(real(z)) && return im*tan(imag(z))
	tanh(float(z))
end

# Arithmetic operators

Base.:(+)(p1::PiTimes,p2::PiTimes) = PiTimes(p1.x + p2.x)
Base.:(+)(p::PiTimes,::Irrational{:π}) = PiTimes(p.x + one(p.x))
Base.:(+)(::Irrational{:π},p::PiTimes) = PiTimes(p.x + one(p.x))

Base.:(-)(p::PiTimes) = PiTimes(-p.x)
Base.:(-)(p1::PiTimes,p2::PiTimes) = PiTimes(p1.x-p2.x)
Base.:(-)(p::PiTimes,::Irrational{:π}) = PiTimes(p.x - one(p.x))
Base.:(-)(::Irrational{:π},p::PiTimes) = PiTimes(one(p.x) - p.x)

Base.:(/)(p1::PiTimes,p2::PiTimes) = p1.x/p2.x
Base.:(/)(p1::PiTimes,y::Real) = PiTimes(p1.x/y)

Base.:(/)(::Irrational{:π},p::PiTimes) = inv(float(p.x))
Base.:(/)(p::PiTimes,::Irrational{:π}) = float(p.x)

Base.:(*)(p1::PiTimes,y::Real) = PiTimes(p1.x*y)
Base.:(*)(y::Real,p1::PiTimes) = PiTimes(p1.x*y)
Base.:(*)(p1::PiTimes,p2::PiTimes) = float(p1) * float(p2)
Base.:(*)(z::Complex{Bool},p::PiTimes) = p*z # switch the orders to get to something that's non-ambiguous

for op in Symbol[:+,:-,:/,:*]
	@eval Base.$op(p::PiTimes,x::AbstractIrrational) = $op(Float64(p),Float64(x))
	@eval Base.$op(x::AbstractIrrational,p::PiTimes) = $op(Float64(x),Float64(p))
end

Base.:(//)(p::PiTimes,n) = PiTimes(p.x//n)

end # module
