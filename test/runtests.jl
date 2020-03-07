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
        for T in (Float16,Float32,Float64,Int8,Int16,Int32,Int64,Int128,BigInt,BigFloat)
            @test one(PiTimes{T}) === true
            @test one(PiTimes{T}(0)) === true
        end
        @test one(Pi) === true
    end
    @testset "zero" begin
        @test zero(Pi) === PiTimes(0) == 0
        @test zero(typeof(Pi)) === PiTimes(0) == 0
        for T in (Float16,Float32,Float64,Int8,Int16,Int32,Int64,Int128)
            @test zero(PiTimes{T}) isa PiTimes{T}
            @test zero(PiTimes{T}(1)) isa PiTimes{T}
            @test zero(PiTimes{T}(1)) === PiTimes{T}(0)
        end
        for T in (BigInt,BigFloat)
            @test zero(PiTimes{T}) isa PiTimes{T}
            @test zero(PiTimes{T}(1)) isa PiTimes{T}
            @test zero(PiTimes{T}(1)) == PiTimes{T}(0)
        end
    end
end

@testset "conversion" begin
    @testset "promote_rule" begin
        for t in (Int8, Int16, Int32, Int64, Int128, Bool, UInt8, UInt16, UInt32, UInt64, UInt128)
                @test promote_rule(PiTimes{Float16},PiTimes{t}) === PiTimes{Float16}
                @test promote_rule(PiTimes{t},PiTimes{Float16}) === PiTimes{Float16}
        end
        for t1 in (Float32, Float64)
            for t2 in (Int8, Int16, Int32, Int64, Bool, UInt8, UInt16, UInt32, UInt64)
                @test promote_rule(PiTimes{t1},PiTimes{t2}) === PiTimes{t1}
                @test promote_rule(PiTimes{t2},PiTimes{t1}) === PiTimes{t1}
            end
        end
        @test promote_rule(PiTimes{Int},Irrational{:π}) === PiTimes{Int}
        @test promote_rule(PiTimes{Float64},Irrational{:π}) === PiTimes{Float64}
        @test promote_rule(Irrational{:π},PiTimes{Float64}) === PiTimes{Float64}
    end
    @testset "convert" begin
        @testset "to float" begin
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
        @testset "from real" begin
            @test convert(PiTimes{Int},π) === PiTimes(1)
            @test convert(PiTimes{Float64},π) === PiTimes(1.0)
            @test convert(PiTimes,π) === PiTimes(1.0)

            @test convert(PiTimes,2) === PiTimes(2/π)
            for t in (Float16,Float32,Float64)
                @test convert(PiTimes{t},2) === PiTimes{t}(2/π)
            end
            convert(PiTimes{BigFloat},2) === PiTimes{BigFloat}(2/π)
        end
        @testset "PiTimes" begin
            for t in (Int8, Int16, Int32, Int64, Int128, Bool, UInt8, UInt16, UInt32, UInt64, UInt128)
                @test convert(PiTimes{Float16},PiTimes{t}(0)) === PiTimes(Float16(0))
            end
            for t1 in (Float32, Float64)
                for t2 in (Int8, Int16, Int32, Int64, Bool, UInt8, UInt16, UInt32, UInt64)
                    @test convert(PiTimes{t1},PiTimes{t2}(0)) === PiTimes(t1(0))
                end
            end
        end
        @testset "Complex" begin
            @test Complex(Pi,Pi) == Pi + im*Pi
            @test Complex(Pi,π) == Pi + im*Pi
            @test Complex(π,Pi) == Pi + im*Pi
        end
    end
end

@testset "arithmetic" begin
    p = PiTimes(3)
    q = PiTimes(4)

    @testset "negation" begin
        @test -p isa PiTimes{Int}
        @test -p === PiTimes(-p.x)
    end 

    @testset "addition" begin
        @test p + p === PiTimes(2p.x)
        @test p + q === PiTimes(p.x + q.x)
    end

    @testset "subtraction" begin
        @test p - p === PiTimes(0)
        @test p - q === PiTimes(p.x - q.x)
    end

    @testset "multiplication" begin
        @testset "PiTimes" begin
            @test p*p ≈ (p.x*π)^2
            @test p*q ≈ (p.x*q.x)*π^2
        end
        @testset "Real" begin
            @testset "Int" begin
                @test p*2 === PiTimes(p.x*2)
                @test 2*p === PiTimes(p.x*2)
                @test 2p === PiTimes(p.x*2)
            end
            @testset "Bool" begin
                @test true*p === p
                @test p*true === p
                @test false*p === zero(p)
                @test p*false === zero(p)
            end
        end
        
        @testset "Complex" begin
            @testset "Complex{Bool}" begin
                @test p * im === PiTimes(0) + PiTimes(p.x)*im 
                @test im* p === PiTimes(0) + PiTimes(p.x)*im
                @test Complex(true,false)*Pi === Pi + im*zero(Pi)
                @test Complex(false,true)*Pi === zero(Pi) + im*Pi
                @test Pi*Complex(false,true) === zero(Pi) + im*Pi
            end
            @testset "Complex{Int}" begin
                @test (1+im)*p === PiTimes(p.x) + PiTimes(p.x)*im
                @test p*(1+im) === PiTimes(p.x) + PiTimes(p.x)*im
            end
        end
    end

    @testset "division" begin
        @testset "PiTimes" begin
            @test p/p === p.x/p.x
            @test p/q === p.x/q.x == 3/4 # works as 1/4 can be stored exactly
            @test q/p === q.x/p.x
        end
        @testset "Real" begin
            @testset "Int" begin
                @test p/2 === PiTimes(p.x/2)
                @test PiTimes(3)/4 === PiTimes(3/4)
            end
            @testset "Bool" begin
                @test p/true === PiTimes(p.x/1)
                @test p/false === PiTimes(Inf)
                @test true/Pi === true/π
                @test false/Pi === 0.0
            end
        end
        @testset "pi" begin
            @test p/π === float(p.x)
            @test PiTimes(7.3)/π === 7.3
            @test π/PiTimes(1/7.3) ≈ 7.3 # within rounding errors
            @test π/PiTimes(1//7) === float(7)
        end
        @testset "Complex" begin
            @testset "Complex{Bool}" begin
                @test p/im === Complex(PiTimes(zero(Float64)),PiTimes(-p.x/1))
                @test p/Complex(true,false) === Complex(PiTimes(p.x/1),PiTimes(-zero(Float64)))
                @test im/Pi === im/π
                @test Complex(true,false)/Pi === Complex(true,false)/π
            end
            @testset "Complex{Int}" begin
                @test p/(1+im) === Complex(PiTimes(1.5),-PiTimes(1.5))
                @test (1+im)/Pi === (1+im)/π
            end
        end
    end

    @testset "Rational" begin
        @test p//2 === PiTimes(p.x//2)
        @test Pi//2 === PiTimes(1//2)
    end

    @testset "Irrational" begin
        @testset "pi" begin
            @test Pi + π === 2Pi
            @test π + Pi === 2Pi
            @test Pi - π === zero(Pi) == 0
            @test π - Pi === zero(Pi) == 0
            @test π*Pi === π^2
            @test Pi*π === π^2
            @test Pi/π === π/π
            @test π/Pi === π/π
            @test Pi^π === π^π
            @test π^Pi === π^π
        end
        @testset "e" begin
            @test Pi + ℯ === π + ℯ
            @test ℯ + Pi === ℯ + π
            @test Pi - ℯ === π - ℯ
            @test ℯ - Pi === ℯ - π
            @test ℯ*Pi === ℯ*π
            @test Pi*ℯ === π*ℯ
            @test Pi/ℯ === π/ℯ
            @test ℯ/Pi === ℯ/π
            @test Pi^ℯ === π^ℯ
            @test ℯ^Pi === ℯ^π
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

@testset "show" begin
    io = IOBuffer()
    show(io,Pi)
    @test String(take!(io)) == "Pi"
    show(io,2Pi)
    @test String(take!(io)) == "2Pi"
end