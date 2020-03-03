module MultiplesOfPi

export PiTimes

struct PiTimes{T<:Real} <: AbstractIrrational
	x :: T
end

Base.:(<)(p::PiTimes,q::PiTimes) = p.x < q.x

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

end # module
