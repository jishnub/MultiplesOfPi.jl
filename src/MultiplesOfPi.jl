module MultiplesOfPi

export PiTimes
export Pi

"""
	PiTimes(x)

Construct a number that behaves as `π*x` for a real `x`
"""
struct PiTimes{T<:Real} <: AbstractIrrational
	x :: T
end

"""
	Pi

The number `PiTimes(1)`, numerically equivalent to `pi`
"""
const Pi = PiTimes(1)

Base.:(<)(p::PiTimes,q::PiTimes) = p.x < q.x
Base.:(==)(p::PiTimes,q::PiTimes) = p.x == q.x
Base.:(==)(p::PiTimes,y::Real) = float(p) == y
Base.:(==)(y::Real,p::PiTimes) = float(p) == y

Base.:(==)(p::PiTimes,::Irrational{:π}) = isone(p.x)
Base.:(==)(::Irrational{:π},p::PiTimes) = isone(p.x)

for T in (:BigFloat,:Float64,:Float32,:Float16)
	@eval Base.$T(p::PiTimes) = π*$T(p.x)
end

for f in (:iszero,:isfinite,:isnan)
	@eval Base.$f(p::PiTimes) = $f(p.x)
end

# Unfortunately Irrational numbers do not have a multiplicative identity of the same type,
# so we make do with something that works
Base.one(::PiTimes{T}) where {T} = one(PiTimes{T})
Base.one(::Type{PiTimes{T}}) where {T} = true

Base.zero(::PiTimes{T}) where {T} = zero(PiTimes{T})
Base.zero(::Type{PiTimes{T}}) where {T} = PiTimes(zero(T))

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
Base.:(*)(z::Complex{Bool},p::PiTimes) = Complex(real(z)*p,imag(z)*p)
Base.:(*)(p::PiTimes,z::Complex{Bool}) = Complex(real(z)*p,imag(z)*p)
Base.:(*)(x::Bool,p::PiTimes) = ifelse(x,p,zero(p))
Base.:(*)(p::PiTimes,x::Bool) = ifelse(x,p,zero(p))

for op in Symbol[:+,:-,:/,:*]
	@eval Base.$op(p::PiTimes,x::AbstractIrrational) = $op(Float64(p),Float64(x))
	@eval Base.$op(x::AbstractIrrational,p::PiTimes) = $op(Float64(x),Float64(p))
end

Base.:(//)(p::PiTimes,n) = PiTimes(p.x//n)
Base.:(//)(n,p::PiTimes) = PiTimes(n//p.x)
Base.:(//)(p::PiTimes,q::PiTimes) = (p.x-q.x)*PiTimes(1)

# Conversion and promotion

for t in (Int8, Int16, Int32, Int64, Int128, Bool, UInt8, UInt16, UInt32, UInt64, UInt128)
    @eval Base.promote_rule(::Type{PiTimes{Float16}}, ::Type{PiTimes{$t}}) = PiTimes{Float16}
    @eval Base.promote_rule(::Type{PiTimes{$t}},::Type{PiTimes{Float16}}) = PiTimes{Float16}
end

for t1 in (Float32, Float64)
    for t2 in (Int8, Int16, Int32, Int64, Bool, UInt8, UInt16, UInt32, UInt64)
        @eval begin
            Base.promote_rule(::Type{PiTimes{$t1}}, ::Type{PiTimes{$t2}}) = PiTimes{$t1}
            Base.promote_rule(::Type{PiTimes{$t2}}, ::Type{PiTimes{$t1}}) = PiTimes{$t1}
        end
    end
end

Base.promote_rule(::Type{PiTimes{T}}, ::Type{Irrational{:π}}) where {T} = PiTimes{T}
Base.promote_rule(::Type{Irrational{:π}}, ::Type{PiTimes{T}}) where {T} = PiTimes{T}

Base.convert(::Type{PiTimes},x::Real) = PiTimes(x/π)
Base.convert(::Type{PiTimes{T}},x::Real) where {T} = PiTimes{T}(x/π)
Base.convert(::Type{PiTimes{T}},p::PiTimes) where {T} = PiTimes{T}(p.x)

function Base.show(io::IO,p::PiTimes)
	pxstr = isone(p.x) ? "" : string(p.x)
	str = iszero(p) ? string(p.x) : pxstr*"Pi"
	print(io,str)
end

end # module
