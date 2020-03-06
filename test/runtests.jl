using MultiplesOfPi
using Test

@testset "Constructor" begin
    @test PiTimes(1) ≈ π
    @test PiTimes(1.5) ≈ π*1.5
    @test PiTimes(1//2) ≈ π/2

    @test PiTimes(1) === Pi
    @test Pi == π
    @test π == Pi
    @test 2π ≈ 2Pi
end

@testset "Ordering" begin
    @testset "inequality" begin
        @test PiTimes(1) < 4
        @test PiTimes(2) > 4
        @test PiTimes(2) > π
        @test PiTimes(2) < PiTimes(4)
    end
    @testset "equality" begin
        @test PiTimes(2) == PiTimes(2)
        @test Pi*0 == 0
        @test 0 == Pi*0
        @test Pi != 2Pi
        @test Pi != 4
    end
end

@testset "zero and one" begin
    @testset "one" begin
        for T in (Float16,Float32,Float64,Int8,Int16,Int32,Int64,Int128)
            @test one(PiTimes{T}) === one(T)
        end
        for T in (BigInt,BigFloat)
            @test one(PiTimes{T}) == one(T)
        end
        @test one(Pi) == 1
    end
    @testset "zero" begin
        @test zero(Pi) == 0
        @test zero(typeof(Pi)) == 0
        for T in (Float16,Float32,Float64,BigFloat,Int8,Int16,Int32,Int64,Int128,BigInt)
            @test zero(PiTimes{T}) isa PiTimes{T}
        end
    end
end

@testset "conversion" begin
    @test convert(Float64,PiTimes(1)) ≈ Float64(1π)
    @test Float64(PiTimes(1)) ≈ Float64(1)*π
    @test convert(Float64,PiTimes(2)) ≈ Float64(2π)
    @test Float64(PiTimes(2)) ≈ Float64(2)*π
    @test convert(BigFloat,PiTimes(1)) ≈ BigFloat(1)*π
    @test BigFloat(PiTimes(1)) ≈ BigFloat(1)*π
    @test convert(Float32,PiTimes(1)) ≈ Float32(1)*π
    @test Float32(PiTimes(1)) ≈ Float32(1)*π
    @test convert(Float16,PiTimes(1)) ≈ Float16(1)*π
    @test Float16(PiTimes(1)) ≈ Float16(1)*π
end

@testset "arithmetic" begin
    p = PiTimes(3)
    q = PiTimes(4)

    @testset "negation" begin
        @test -p isa PiTimes{Int}
        @test -p == PiTimes(-p.x)
    end 

    @testset "addition" begin
        @test p + p == PiTimes(2p.x)
        @test p + q == PiTimes(p.x + q.x)
    end

    @testset "subtraction" begin
        @test p - p == PiTimes(0)
        @test p - q == PiTimes(p.x - q.x)
    end

    @testset "multiplication" begin
        @test p*2 == PiTimes(p.x*2)
        @test 2*p == PiTimes(p.x*2)
        @test 2p == PiTimes(p.x*2)
        @test p*p ≈ (p.x*π)^2
        @test p*q ≈ (p.x*q.x)*π^2
        @test p * im == PiTimes(0) + PiTimes(p.x)*im 
        @test im* p == PiTimes(0) + PiTimes(p.x)*im
        @test (1+im)*p == PiTimes(p.x) + PiTimes(p.x)*im
    end

    @testset "division" begin
        @test p/p == 1
        @test p/q == p.x/q.x == 3/4 # works as 1/4 can be stored exactly
        @test q/p == q.x/p.x

        @test p/2 == PiTimes(p.x/2)

        @test p/π === float(p.x)
        @test PiTimes(7.3)/π === 7.3
        @test π/PiTimes(1/7.3) ≈ 7.3 # within rounding errors
        @test π/PiTimes(1//7) === float(7)

        @test PiTimes(3)/4 == PiTimes(3/4)
    end

    @testset "Rational" begin
        @test p//2 == PiTimes(p.x//2)
        @test Pi//2 == PiTimes(1//2)
    end

    @testset "irrational" begin
        @testset "pi" begin
            @test Pi + π == 2Pi
            @test π + Pi == 2Pi
            @test Pi - π == 0
            @test π - Pi == 0
            @test π*Pi == π^2
            @test Pi*π == π^2
            @test Pi/π == 1
            @test π/Pi == 1
            @test Pi^π == π^π
            @test π^Pi == π^π
        end
        @testset "e" begin
            @test Pi + ℯ == π + ℯ
            @test ℯ + Pi == ℯ + π
            @test Pi - ℯ == π - ℯ
            @test ℯ - Pi == ℯ - π
            @test ℯ*Pi == ℯ*π
            @test Pi*ℯ == π*ℯ
            @test Pi/ℯ == π/ℯ
            @test ℯ/Pi == ℯ/π
            @test Pi^ℯ == π^ℯ
            @test ℯ^Pi == ℯ^π
        end
    end
end

@testset "trigonometric functions" begin
    @testset "sin and cos" begin
        x = 1.5
        @test sin(PiTimes(x)) == sinpi(x) == -1 ≈ sin(π*x)
        x = 1
        @test cos(PiTimes(x)) == cospi(x) == -1 ≈ cos(π*x)
        
        a = sincos(PiTimes(x))
        @test a[1] == sinpi(x)
        @test a[2] == cospi(x)
    end

    @testset "complex exp" begin
        x = 1.1
        @test cis(PiTimes(x)) == cospi(x) + im*sinpi(x) ≈ cos(π*x) + im*sin(π*x)

        @test exp(im*PiTimes(0)) == cis(PiTimes(0)) == 1
        @test exp(im*PiTimes(1/2)) == cis(PiTimes(1/2)) == im
        @test exp(im*PiTimes(1)) == cis(PiTimes(1))== -1
        @test exp(im*PiTimes(3/2)) == cis(PiTimes(3/2)) == -im

        @test exp(im*Pi/2) == im
        @test exp(im*Pi) == -1
        @test exp(im*3Pi/2) == -im
    end

    @testset "sec and csc" begin
        @test sec(PiTimes(1)) == -1
        @test sec(Pi) == -1
        @test sec(PiTimes(2)) == 1
        @test sec(2Pi) == 1
        @test !isfinite(sec(PiTimes(1/2)))
        @test !isfinite(sec(Pi/2))

        @test csc(PiTimes(1/2)) == 1
        @test csc(Pi/2) == 1
        @test csc(PiTimes(3/2)) == -1
        @test csc(3Pi/2) == -1
        @test !isfinite(csc(PiTimes(1))) 
        @test !isfinite(csc(Pi))
    end

    @testset "tan and cot" begin
        @testset "tan" begin
            @test isinf(tan(Pi/2))
            @test isinf(tan(-Pi/2))
            @test tan(Pi/4) == 1
            @test tan(-Pi/4) == -1
            @test tan(Pi*0) == 0
            @test tan(-Pi*0) == 0
            @test tan(Pi) == 0
            @test tan(-Pi) == 0
            @test tan(2Pi) == 0
            @test tan(-2Pi) == 0
            @test tan(Pi/8) == tan(pi/8)

            @test tan(nextfloat(0.0)*π) == tan(PiTimes(nextfloat(0.0)))
            @test tan(prevfloat(0.0)*π) == tan(PiTimes(prevfloat(0.0)))
        end
        @testset "cot" begin
            @test cot(Pi/2) == 0
            @test cot(-Pi/2) == 0
            @test cot(Pi/4) == 1
            @test cot(-Pi/4) == -1
            @test isinf(cot(Pi*0))
            @test isinf(cot(-Pi*0))
            @test isinf(cot(Pi))
            @test isinf(cot(-Pi))
            @test isinf(cot(2Pi))
            @test isinf(cot(-2Pi))
            @test cot(Pi/8) == cot(pi/8)

            @test cot(nextfloat(0.0)*π) == cot(PiTimes(nextfloat(0.0)))
            @test cot(prevfloat(0.0)*π) == cot(PiTimes(prevfloat(0.0)))
        end
    end

    @testset "sinc" begin
        @test sinc(Pi) == sinc(π)
    end
end

@testset "hyperbolic" begin
    @testset "cosh" begin
        @test cosh(im*PiTimes(1/2)) == 0
    end
    @testset "sinh" begin
        @test sinh(im*PiTimes(-1)) == 0
    end
    @testset "tanh" begin
        @test tanh(im*Pi/4) == im
        @test tanh(im*Pi/2) == im*Inf
        @test tanh(im*Pi/3) == im*tan(Pi/3)
        @test tanh((1+im)*Pi) == tanh((1+im)*π)
        @test tanh(2+im*Pi) == tanh(2+im*π)
    end
end

@testset "Inf NaN zero check" begin
    @test isinf(PiTimes(Inf))
    @test !isfinite(PiTimes(Inf))
    @test isfinite(Pi)
    @test isnan(PiTimes(NaN))
    @test iszero(PiTimes(0))
end