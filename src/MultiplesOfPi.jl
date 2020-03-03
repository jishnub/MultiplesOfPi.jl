module MultiplesOfPi

export PiTimes
export Pi

struct PiTimes{T<:Real} <: AbstractIrrational
	x :: T
end

const Pi = PiTimes(1)

Base.:(<)(p::PiTimes,q::PiTimes) = p.x < q.x
Base.:(==)(p::PiTimes,q::PiTimes) = p.x == q.x

Base.:(==)(p::PiTimes,::Irrational{:π}) = isone(p.x)
Base.:(==)(::Irrational{:π},p::PiTimes) = isone(p.x)

for T in (:BigFloat,:Float64,:Float32,:Float16)
	eval(quote
		Base.$T(p::PiTimes) = π*$T(p.x)
	end)
end

# Define trigonometric functions

@inline Base.sin(p::PiTimes) = sinpi(p.x)
@inline Base.cos(p::PiTimes) = cospi(p.x)
@inline Base.sincos(p::PiTimes) = (sin(p),cos(p))

@inline Base.cis(p::PiTimes) = Complex(cos(p),sin(p))

Base.:(+)(p1::PiTimes,p2::PiTimes) = PiTimes(p1.x+p2.x)

Base.:(-)(p::PiTimes) = PiTimes(-p.x)
Base.:(-)(p1::PiTimes,p2::PiTimes) = PiTimes(p1.x-p2.x)

Base.:(/)(p1::PiTimes,p2::PiTimes) = p1.x/p2.x
Base.:(/)(p1::PiTimes,y::Real) = PiTimes(p1.x/y)

Base.:(/)(::Irrational{:π},p::PiTimes) = inv(float(p.x))
Base.:(/)(p::PiTimes,::Irrational{:π}) = float(p.x)

Base.:(*)(p1::PiTimes,y::Real) = PiTimes(p1.x*y)
Base.:(*)(y::Real,p1::PiTimes) = PiTimes(p1.x*y)
Base.:(*)(p1::PiTimes,p2::PiTimes) = float(p1) * float(p2)
Base.:(*)(z::Complex{Bool},p::PiTimes) = p*z # switch the orders to get to something that's non-ambiguous

Base.:(//)(p::PiTimes,n) = PiTimes(p.x//n)

Base.one(::Type{PiTimes{T}}) where {T} = one(T)

end # module
