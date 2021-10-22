using MultiplesOfPi
using Test

@testset "Constructor" begin
    @test (@inferred (1)*Pi) ≈ π
    @test (@inferred (1.5)*Pi) ≈ 1.5*pi
    @test (@inferred (1//2)*Pi) ≈ π/2

    @test (1)*Pi === Pi
    @test Pi == π
    @test π == Pi
    @test 2π ≈ 2Pi

    @test float(Pi^0.5) == pi^0.5

    @testset "nested" begin
        @test PiExpTimes{PiExpTimes{Int}}(3, 0) == 3
        @test PiExpTimes{PiExpTimes{Float64}}(3, 0) == PiExpTimes{Float64}(3)
        @test PiExpTimes{PiExpTimes{Float64}}(Pi, 0) == Pi

        @test PiExpTimes{PiExpTimes}(Pi) isa PiExpTimes{PiExpTimes}
        @test PiExpTimes{PiExpTimes}(Pi) == Pi
        @test PiExpTimes{PiExpTimes{Int}}(Pi) isa PiExpTimes{PiExpTimes{Int}}
        @test PiExpTimes{PiExpTimes{Int}}(Pi) == Pi
    end

    @testset "Irrational{:π}" begin
        @test (π)*Pi == PiExpTimes(1, 2)
        @test PiExpTimes{Irrational{:π}}(π) isa PiExpTimes{Irrational{:π}}
        @test PiExpTimes{Irrational{:π}}(π).x === π

        @test PiExpTimes{Real}(π) isa PiExpTimes{Real}
        @test PiExpTimes{Real}(π).x === 1

        @test PiExpTimes{Irrational{:π}}(π, 2) == PiExpTimes(1, 3)
    end

    @test PiExpTimes(1, -1) ≈ 1/π
    @test PiExpTimes(2, -2) ≈ 2/π^2

    @test PiExpTimes(Pi^2, 2) === Pi^4
    @test PiExpTimes(Pi^2, -2) == 1
    @test PiExpTimes(1.0Pi^2, -2) == 1.0
    @test PiExpTimes(Pi^-2, 2) == 1.0
    @test PiExpTimes{Int}(Pi^2, 2) === Pi^4
    @test PiExpTimes{Float64}(Pi^2, 2) === PiExpTimes(1.0, 4)

    @test PiExpTimes{Int}(2.0, 2) === PiExpTimes(2, 2)

    @test (pi)*Pi === PiExpTimes(1, 2)
    @test PiExpTimes(pi, 2) === PiExpTimes(1, 3)

    @test PiExpTimes(π, 2) === PiExpTimes(1, 3)

    @test PiExpTimes(PiExpTimes(1, 3), 2) === Pi^5

    @test float(PiExpTimes(2, 2)) === 2*pi*pi

    @test float(PiExpTimes(BigFloat(2)^10000)) isa BigFloat
end

@testset "Ordering" begin
    @testset "inequality" begin
        @testset "PiExpTimes" begin
            @test PiExpTimes(3, 2) < PiExpTimes(4, 2)
            @test PiExpTimes(3, 1) < PiExpTimes(3, 2)
            @test PiExpTimes(3, 1) < PiExpTimes(3.5, 2)
            @test (3)*Pi < PiExpTimes(3, 2)
            @test π < PiExpTimes(1, 2)
            @test !(PiExpTimes(1, 2) < π)
            @test 1 < PiExpTimes(1, 2)
            @test PiExpTimes(1, 2) > π
            @test PiExpTimes(1, 2) > 1
            @test PiExpTimes(1, 0) < 2
            @test PiExpTimes(1, 0) < π
        end
    end
    @testset "equality" begin
        @testset "PiTimes" begin
            @test Pi == π
            @test π == Pi
            @test (2)*Pi == (2)*Pi
            @test (0)*Pi == 0
            @test 0 == (0)*Pi
            @test Pi != 2Pi
            @test Pi != 4
        end
        @testset "PiExpTimes" begin
            Pi² = PiExpTimes(1, 2)
            Pi³ = PiExpTimes(1, 3)

            z² = PiExpTimes(0, 2)
            z³ = PiExpTimes(0, 3)

            @test Pi² == Pi²

            @test Pi² != π^2
            @test π^2 != Pi²

            @test Pi² != 3
            @test 3 != Pi²

            @test Pi² != Pi
            @test Pi  != Pi²

            @test Pi² != Pi³
            @test Pi³ != Pi²

            @test !(Pi² == π)
            @test !(π == Pi²)

            @test z² == z³
            @test z² == 0
        end
        @testset "nested" begin
            @test π == PiExpTimes{PiExpTimes{Int}}(Pi^2, -1)
            @test PiExpTimes{PiExpTimes{Int}}(Pi^2, -1) == π
        end
    end
end

@testset "sign" begin
    @testset "PiExpTimes" begin
        @testset "sign" begin
            @test sign(PiExpTimes(1, 2)) == 1
            @test sign(PiExpTimes(3.4, 2)) == 1.0
            @test sign(PiExpTimes(0.0, 2)) == 0.0
            @test sign(PiExpTimes(-0.0, 2)) == -0.0
            @test sign(PiExpTimes(1//2, 2)) == 1//1
        end
        @testset "signbit" begin
            @test !signbit(PiExpTimes(1, 2))
            @test !signbit(PiExpTimes(0, 2))
            @test !signbit(PiExpTimes(-0, 2))
            @test signbit(PiExpTimes(-1, 2))
            @test signbit(PiExpTimes(-0.0, 2))
            @test !signbit(PiExpTimes(0.0, 2))
        end
    end
end

@testset "zero and one" begin
    @testset "one" begin
        for T in (Float16,Float32,Float64,Int8,Int16,Int32,Int64,Int128,BigInt,BigFloat)
            @test one(PiExpTimes{T,Int}) == 1
            @test one(PiExpTimes{T,Float64}) == 1
            @test one(PiExpTimes{T,T}) == 1
            @test one(PiExpTimes(zero(T))) == 1
        end

        @test one(Pi) * Pi == Pi
        @test one(Pi^2) * Pi^2 == Pi^2

        @test one(PiExpTimes{PiExpTimes{Int}}(1)) == 1
    end
    @testset "zero" begin
        @test zero(Pi) == 0
        @test zero(typeof(Pi)) == 0

        @test zero(Pi) + Pi == Pi
        @test zero(Pi^2) + Pi^2 == Pi^2

        for T in (Float16,Float32,Float64,Int8,Int16,Int32,Int64,Int128)
            @test zero(PiExpTimes{T,Int}) isa PiExpTimes{T,Int}
            @test zero(PiExpTimes{T,Float64}) isa PiExpTimes{T,Float64}
            @test zero(PiExpTimes{T,T}) isa PiExpTimes{T,T}
            @test zero(PiExpTimes{T}(1)) isa PiExpTimes{T}
            @test zero(PiExpTimes{T}(1)) == 0
        end
        for T in (BigInt,BigFloat)
            @test zero(PiExpTimes{T,Int}) isa PiExpTimes{T,Int}
            @test zero(PiExpTimes{T,BigInt}) isa PiExpTimes{T,BigInt}
            @test zero(PiExpTimes{T}(1)) isa PiExpTimes{T}
            @test zero(PiExpTimes{T}(1)) == 0
        end
    end
end

@testset "Range" begin
    @testset "LinRange" begin
        l = LinRange(0.0Pi,1.0Pi,11)
        lpi = LinRange(0,π,11)
        @test l[1] === 0.0Pi
        @test l[end] === 1.0Pi

        for i in eachindex(l,lpi)
            @test l[i] ≈ lpi[i]
        end
    end
    @testset "range" begin
        @testset "step" begin
            r = range(0Pi, step=Pi, length = 3)
            @test length(r) == 3
            @test r[1] == 0Pi
            @test r[end] == 2Pi
            @test step(r) == 1Pi
        end
    end
end

@testset "arithmetic" begin
    p = 3Pi
    q = 4Pi
    r = 4Pi^2

    @testset "negation" begin
        @test -p isa PiExpTimes{Int}
        @test -r isa PiExpTimes{Int}
        @test -p === -p.x*Pi
        @test -r === PiExpTimes(-r.x, 2)
    end

    @testset "addition" begin
        @test p + p === (2p.x)*Pi
        @test p + q === (p.x + q.x)*Pi
        @test r + r === PiExpTimes(2r.x, 2)
        @test p + r === (p.x + r.x*π)*Pi
        @test p + r ≈ float(p) + float(r)
    end

    @testset "subtraction" begin
        @test p - p === (0)*Pi == 0
        @test p - q === (p.x - q.x)*Pi
        @test r - r === PiExpTimes(zero(r.x), 2) == 0
        @test p - r === (p.x - r.x*π)*Pi
        @test p - r ≈ float(p) - float(r)
    end

    @testset "multiplication" begin
        @testset "PiTimes" begin
            @test (@inferred p*p) === PiExpTimes(p.x^2, 2)
            @test (@inferred p*p) ≈ (p.x*π)^2
            @test (@inferred p*q) === PiExpTimes(p.x*q.x, 2)
            @test (@inferred q*p) === PiExpTimes(p.x*q.x, 2)
            @test (@inferred p*q) ≈ (p.x*q.x)*π^2
            @test (@inferred p*r) === PiExpTimes(p.x*r.x, 3) ≈ (p.x*r.x)*π^3
            @test (@inferred r*p) === PiExpTimes(p.x*r.x, 3) ≈ (p.x*r.x)*π^3
        end
        @testset "Real" begin
            @testset "Int" begin
                @test p*2 === (p.x*2)*Pi
                @test 2*p === (p.x*2)*Pi
                @test 2p === (p.x*2)*Pi
                @test 2r === PiExpTimes(2r.x, 2)

                @test PiExpTimes(3, 0) * 2 == 6
                @test 2 * PiExpTimes(3, 0) == 6
            end
            @testset "Bool" begin
                @test (@inferred true*p) === p
                @test (@inferred p*true) === p
                @test (@inferred false*p) === zero(p)
                @test (@inferred p*false) === zero(p)

                @test PiExpTimes(3, 0) * true == 3
                @test PiExpTimes(3, 0) * false == 0
                @test PiExpTimes(3.0, 0) * true == 3.0
                @test PiExpTimes(3.0, 0) * false == 0.0
                @test PiExpTimes(-3.0, 0) * true == -3.0
                @test PiExpTimes(-3.0, 0) * false == -0.0
                @test sign(PiExpTimes(-3.0, 0) * false) === -0.0
                @test -3.0Pi * false === false * (-3.0Pi)
            end
            @testset "float" begin
                @test (@inferred p*1/p) == 1.0
                @test (@inferred r*1/r) == 1.0

                @test 2.3*p === (2.3*p.x)*Pi

                @test PiExpTimes(3.0, 0) * 2 == 6.0
                @test 2 * PiExpTimes(3.0, 0) == 6.0

                @test 2 * PiExpTimes(-3.0, 0) == -6.0
            end
        end

        @testset "Complex" begin
            @testset "Complex{Bool}" begin
                @test (@inferred p * im) === (0)*Pi + (p.x)*Pi*im
                @test (@inferred im* p) === (0)*Pi + (p.x)*Pi*im
                @test Complex(true,false)*Pi === Pi + im*zero(Pi)
                @test Complex(false,true)*Pi === zero(Pi) + im*Pi
                @test Pi*Complex(false,true) === zero(Pi) + im*Pi
            end
            @testset "Complex{Int}" begin
                @test (1+im)*p === (p.x)*Pi + (p.x)*Pi*im
                @test p*(1+im) === (p.x)*Pi + (p.x)*Pi*im
            end
        end
    end

    @testset "division" begin
        @testset "PiTimes" begin
            @test p/p == p.x/p.x
            @test p/q == p.x/q.x == 3/4 # works as 1/4 can be stored exactly
            @test q/p == q.x/p.x
        end
        @testset "PiExpTimes" begin
            @test r/r == 1.0
            @test r/p == (r.x/p.x)*Pi
            @test p/r == PiExpTimes(p.x/r.x, -1)
        end
        @testset "Real" begin
            @testset "Int" begin
                @test p/2 === (p.x/2)*Pi
                @test (3)*Pi/4 === (3/4)*Pi
            end
            @testset "Bool" begin
                @test p/true == (p.x/true)*Pi
                @test p/false == (Inf)*Pi
                @test true/Pi == PiExpTimes(1.0, -1)
                @test false/Pi == PiExpTimes(0.0, -1)
            end
        end
        @testset "pi" begin
            @test p/π == p.x
            @test (7.3)*Pi/π == 7.3
            @test π/((1/7.3)*Pi) ≈ 7.3 # within rounding errors
            @test π/((1//7)*Pi) == 7//1

            @test Pi^2/π == Pi
            @test π/Pi^2 == Pi^-1
        end
        @testset "Complex" begin
            @testset "Complex{Bool}" begin
                @test p/im === Complex(zero(Float64)*Pi,(-p.x/1)*Pi)
                @test p/Complex(true,false) == Complex((p.x/1),-zero(Float64))*Pi
                @test im/Pi === im*PiExpTimes(1.0, -1)
                @test Complex(true,false)/Pi == Complex(true,false)*PiExpTimes(1.0, -1)
                @test Pi*Pi/im == -1.0Pi*Pi*im
            end
            @testset "Complex{Int}" begin
                @test p/(1+im) === Complex(1.5Pi,-1.5Pi)
                @test Pi/2im === Complex(0.0Pi,-0.5Pi)
                @test (1+im)/Pi === (1+im)*PiExpTimes(1.0, -1)
            end
            @testset "Complex{PiTimes}" begin
                @testset "Real/Complex{PiExpTimes}" begin
                    @test Pi*Pi/(im*Pi*Pi) == Complex(0.0,-1.0)
                    resim = -PiExpTimes(1.0, -1)
                    @test Pi/(im*Pi*Pi) == Complex(zero(resim),resim)
                end

                @testset "Complex{PiExpTimes}/Real" begin
                    @testset "Complex{PiExpTimes}/Int" begin
                        @test Pi*Pi*im/1 == Complex(0.0*Pi^2,1.0*Pi^2)
                    end
                    @testset "Complex{PiExpTimes}/PiExpTimes" begin
                        @test (Pi*im)/Pi == 1.0im
                        @test Pi*Pi*im/Pi == 1.0Pi*im
                        @test Pi*Pi*im/(Pi*Pi) == Complex(0.0,1.0)
                    end
                end

                @testset "Complex{PiExpTimes}/Complex" begin
                    @testset "Complex{PiExpTimes}/Complex{Int}" begin
                        z = Pi*(1+2im)/(1+2im)
                        @test real(z) == (1.0)*Pi
                        @test iszero(imag(z))
                    end
                    @testset "Complex{PiExpTimes}/Complex{PiExpTimes}" begin
                        z = Pi*(1+2im)
                        @test isone(real(z/z))
                        @test iszero(imag(z/z))

                        @test Pi^2*im/(im*Pi^2) == Complex(1.0,0.0)
                        @test Pi^2*im/(im*Pi) == Complex(1.0,0.0)*Pi
                        @test Pi*im/(im*Pi^2) == Complex(1.0,0.0)*Pi^-1
                    end
                end
            end
        end
    end

    @testset "exponentiation" begin
        @testset "Int exponent" begin
            @test Pi^0 === PiExpTimes(1, 0) == 1
            @test Pi^1 === Pi == π
            @test Pi^2 === PiExpTimes(1, 2) ≈ π*π
            @test Pi^-1 === 1/Pi === PiExpTimes(1.0, -1) ≈ 1/π
            @test Pi^-2 === 1/Pi^2 === PiExpTimes(1.0, -2) ≈ 1/π^2

            for n=1:10
                @test Pi^n ≈ π^n
                @test Pi^-n == 1/Pi^n ≈ float(π)^-n
                @test float(Pi^-n) === float(1/Pi^n) ≈ 1/π^n
            end

            @test (Pi^2)^0 == 1

            @test float(Pi^big(10000)) == BigFloat(pi)^10000
        end

        @testset "Float64 exponent" begin
            @test float(@inferred Pi^2.0) === pi^2.0
            @test float(@inferred Pi^-2.0) === pi^-2.0
            @test (@inferred Pi^0.0) == 1.0
            @test float(@inferred Pi^0.5) == float(@inferred sqrt(Pi))
            @test BigFloat(@inferred Pi^0.5) == BigFloat(@inferred sqrt(Pi))
            @test float(@inferred Pi^(1//3)) == float(@inferred cbrt(Pi))
            @test BigFloat(@inferred Pi^(1//3)) == BigFloat(@inferred cbrt(Pi))
            for T in [Float64, BigFloat]
                @test T(sqrt(Pi)^2) == T(Pi)
                @test T(sqrt(Pi^2)) == T(Pi)
                @test T(sqrt(Pi^2.0)) == T(Pi)
                @test T(sqrt(Pi^(2//1))) == T(Pi)
                @test T(cbrt(Pi)^3) == T(Pi)
                @test T(cbrt(Pi^3)) == T(Pi)
                @test T(cbrt(Pi^3.0)) == T(Pi)
                @test T(cbrt(Pi^(3//1))) == T(Pi)
            end
        end

        @testset "Irrational exponent" begin
            for T in [Float64, BigFloat]
                @test T(@inferred Pi^Pi) == T(pi)^pi
                @test T(@inferred Pi^π) == T(pi)^pi
                @test T(@inferred π^Pi) == T(pi)^pi
            end
        end

        @testset "Rational exponent" begin
            @test Float64(@inferred Pi^(3//4)) == pi^(3//4)
        end

        @testset "Pi as exponent" begin
            @test (@inferred 1^Pi) == 1
            @test (@inferred 2^Pi) == 2^π
            @test (@inferred 2.0^Pi) == 2.0^π
            @test (@inferred 2.0^(1.0Pi)) == 2.0^(1.0*π)
            @test (@inferred big(2)^Pi) == big(2)^π
            @test (@inferred big(2.0)^Pi) == big(2.0)^π
        end

        @testset "reciprocal" begin
            @test (@inferred Pi * Pi^-1)== 1.0
            @test (@inferred Pi^-1 * Pi)== 1.0
            @test (@inferred Pi^-1.0 * Pi)== 1.0
            @test (@inferred 2Pi * Pi^-1)== 2.0
            @test (@inferred 2Pi^2 * Pi^-2) == 2.0
            @test (@inferred 2Pi^2.0 * Pi^-2) == 2.0
        end
    end

    @testset "inv" begin
        @test inv(Pi) === PiExpTimes(1.0, -1)
        @test inv(Pi^2) === PiExpTimes(1.0, -2)
        @test inv(PiExpTimes{Int}(1, 0)) == 1.0
    end

    @testset "Rational" begin
        @test p//2 === (p.x//2)*Pi
        @test Pi//2 === (1//2)*Pi
        @test Pi//Pi == 1//1
        @test 3Pi//4Pi == 3//4
        @test (1//2)Pi//3Pi == 1//6
        @test 2Pi//(3//2)Pi == 4//3
        @test (1//2)Pi//(1//3)Pi == 3//2
        @test (1//2 + 3//2 * im)Pi == (1+3im)*Pi // 2

        @test PiExpTimes{Int}(1, 0)//2 == 1//2

        @testset "pi" begin
            @test π//Pi == 1//1
            @test Pi//π == 1//1

            @test π//Pi^2 == (1//Pi)
            @test Pi^2//π == Pi//1
            @test PiExpTimes{Int}(1, 0)//Pi == 1//Pi
        end

        @testset "PiExpTimes // PiExpTimes" begin
            @test Pi^2//Pi === Pi//1
            @test Pi//Pi^2 === 1//Pi
            @test (3Pi^2//2)//(Pi//2) == 3Pi//1
            @test (3Pi//2)//(Pi^2//2) == (3//Pi)
        end

        @testset "Real // Complex" begin
            @testset "Complex{Int}" begin
                @test Pi//(2im) === (0//1 - 1//2*im)*Pi
                @test Pi^2//(2im) === (0//1 - 1//2*im)*Pi^2
                @test Pi^-1//(2im) === (0//1 - 1//2*im)*Pi^-1
            end
            @testset "Complex{PiTimes}" begin
                @test 2Pi // (2Pi + 0Pi*im) == 1//1 + 0//1*im
                @test 2Pi^2 // (2Pi^2 + 0Pi^2*im) == 1//1 + 0//1*im
                @test 2Pi^2 // (2Pi + 0Pi*im) === (1//1)*Pi + (0//1)*Pi*im
                @test 2Pi // (2Pi^2 + 0Pi^2*im) === (1//Pi) + (0//Pi)*im
            end
        end

        @testset "Complex // Real" begin
            @test (2//1 + 1//1*im)Pi // 2 === (1//1 + 1//2*im)Pi
            @test (2//1 + 1//1*im)Pi^2 // 2 === (1//1 + 1//2*im)Pi^2
            @test (2//1 + 1//1*im)Pi // Pi == (2//1 + 1//1*im)
            @test (2//1 + 1//1*im)Pi^2 // Pi^2 == (2//1 + 1//1*im)
            @test (2//1 + 1//1*im)Pi^2 // Pi === (2//1 + 1//1*im)*Pi
            @test (2//1 + 1//1*im)Pi // Pi^2 === (2//Pi + (1//Pi)*im)
            @test Complex(π,π)//Pi == Complex(1//1,1//1)
            @test Complex(π,π)//Pi^2 === Complex(1//Pi,1//Pi)
        end

        @testset "Complex // Complex" begin
            z1 = (1//2 + 3//2 * im)Pi
            z2 = (1//2 + 3//2 * im)Pi^2
            @test z1 // z1 == (1//1 + 0//1 *im)
            @test z2 // z2 == (1//1 + 0//1 *im)
            @test z2 // z1 === (1//1 + 0//1 * im)*Pi
            @test z1 // z2 === (1 + 0*im)//Pi
        end
    end

    @testset "Irrational" begin
        @testset "pi" begin
            @test Pi + π == 2Pi
            @test Pi^1.0 + π == 2Pi
            @test π + Pi == 2Pi
            @test Pi - π == zero(Pi) == 0
            @test π - Pi == zero(Pi) == 0
            @test π*Pi == Pi*Pi
            @test Pi*π == Pi*Pi
            @test 2Pi*π == (2Pi)*Pi
            @test π*2Pi == (2Pi)*Pi
            @test Pi/π == 1
            @test π/Pi == 1.0
            @test (Pi)*Pi/π == Pi
            @test Pi*Pi/π == Pi
            @test Pi*π/Pi == (one(Float64))*Pi
            @test Pi*π/Pi == Pi == pi
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
            @test float(Pi^ℯ) == pi^ℯ
        end
    end
end

@testset "trigonometric functions" begin
    @testset "PiTimes" begin
        @testset "sin and cos" begin
            x = 1.5
            @test sin(x*Pi) == sinpi(x) == -1 ≈ sin(π*x)
            x = 1
            @test cos(x*Pi) == cospi(x) == -1 ≈ cos(π*x)

            a = sincos(x*Pi)
            @test a[1] == sinpi(x)
            @test a[2] == cospi(x)
        end

        @testset "tan and cot" begin
            @testset "tan" begin
                @test isinf(tan(Pi/2))
                @test isinf(tan(-Pi/2))
                @test (@inferred tan(Pi/4)) == 1
                @test (@inferred tan(-Pi/4)) == -1
                @test (@inferred tan(Pi*0)) == 0
                @test (@inferred tan(-Pi*0)) == 0
                @test (@inferred tan(Pi)) == 0
                @test (@inferred tan(-Pi)) == 0
                @test (@inferred tan(2Pi)) == 0
                @test (@inferred tan(-2Pi)) == 0
                @test (@inferred tan(Pi/8)) == tan(pi/8)

                @test tan(nextfloat(0.0)*π) == tan(nextfloat(0.0)*Pi)
                @test tan(prevfloat(0.0)*π) == tan(prevfloat(0.0)*Pi)
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

                @test cot(nextfloat(0.0)*π) == cot(nextfloat(0.0)*Pi)
                @test cot(prevfloat(0.0)*π) == cot(prevfloat(0.0)*Pi)
            end
        end

        @testset "complex exp" begin
            x = 1.1
            @test cis(x*Pi) == cospi(x) + im*sinpi(x) ≈ cos(π*x) + im*sin(π*x)

            @test exp(im*(0)*Pi) == cis((0)*Pi) == 1
            @test exp(im*(1/2)*Pi) == cis((1/2)*Pi) == im
            @test exp(im*(1)*Pi) == cis((1)*Pi)== -1
            @test exp(im*(3/2)*Pi) == cis((3/2)*Pi) == -im

            @test exp(im*Pi/2) == im
            @test exp(im*Pi) == -1
            @test exp(im*3Pi/2) == -im
        end

        @testset "sec and csc" begin
            @test sec((1)*Pi) == -1
            @test sec(Pi) == -1
            @test sec((2)*Pi) == 1
            @test sec(2Pi) == 1
            @test !isfinite(sec((1/2)*Pi))
            @test !isfinite(sec(Pi/2))

            @test csc((1/2)*Pi) == 1
            @test csc(Pi/2) == 1
            @test csc((3/2)*Pi) == -1
            @test csc(3Pi/2) == -1
            @test !isfinite(csc((1)*Pi))
            @test !isfinite(csc(Pi))
        end

        @testset "sinc" begin
            @test sinc(Pi) === sinc(π)
        end
    end
    @testset "PiExpTimes" begin
        @testset "Pi^2" begin
            @test sin(Pi^2) ≈ sin(π^2)
            @test cos(Pi^2) ≈ cos(π^2)
            s1, c1 = sincos(Pi^2)
            s2, c2 = sincos(π^2)
            @test (s1 ≈ s2) && (c1 ≈ c2)
            @test tan(Pi^2) ≈ tan(π^2)
            @test cot(Pi^2) ≈ cot(π^2)
            @test sinc(Pi^2) ≈ sinc(π^2)
            @test sinc(Pi^-2) ≈ sinc((1/π)^2)
            @test sinc(1/Pi^2) ≈ sinc((1/π)^2)
            @test sec(Pi^2) ≈ sec(π^2)
            @test csc(Pi^2) ≈ csc(π^2)
        end
        @testset "Pi^0" begin
            x = 1
            p = PiExpTimes{Int}(x, 0)
            @test sin(p) == sin(x)
            @test cos(p) == cos(x)
            @test sincos(p) == sincos(x)
            @test tan(p) == tan(x)
            @test cot(p) == cot(x)
            @test sinc(p) == sinc(x)
            @test sinc(p) == sinc(x)
            @test sec(p) == sec(x)
            @test csc(p) == csc(x)
        end
    end
end

@testset "hyperbolic" begin
    @testset "cosh" begin
        @test cosh(im*(1/2)*Pi) == 0
    end
    @testset "sinh" begin
        @test sinh(im*(-1)*Pi) == 0
    end

    @testset "tanh" begin
        @test tanh(im*Pi/4) == im
        @test tanh(im*Pi/2) == im*Inf
        @test tanh(im*Pi/3) == im*tan(Pi/3)
        @test tanh((1+im)*Pi) == tanh((1+im)*π)
        @test tanh(2+im*Pi) == tanh(2+im*π)
    end

    @testset "PiExpTimes" begin
        @testset "Pi^2" begin
            @test sinh(Pi^2) === sinh(π^2)
            @test cosh(Pi^2) === cosh(π^2)
            @test tanh(Pi^2) === tanh(π^2)
        end
        @testset "Pi^0" begin
            x = 1
            p = PiExpTimes{Int}(x, 0)
            @test sinh(p) === sinh(x)
            @test cosh(p) === cosh(x)
            @test tanh(p) === tanh(x)
        end
    end
end

@testset "Inf NaN zero check" begin
    @test isinf((Inf)*Pi)
    @test !isfinite((Inf)*Pi)
    @test isfinite(Pi)
    @test !isfinite(PiExpTimes(Inf, 2))
    @test !isfinite(PiExpTimes(Inf, 0))
    @test !isfinite(PiExpTimes(Inf, -1))
    @test isinf(PiExpTimes(Inf, 2))
    @test isinf(PiExpTimes(Inf, 0))
    @test isinf(PiExpTimes(Inf, -1))
    @test isinf(PiExpTimes(1,Inf))
    @test isfinite(PiExpTimes(1,-Inf))

    @test isnan((NaN)*Pi)
    @test isnan(PiExpTimes(1,NaN))
    @test isnan(PiExpTimes(NaN,1))
    @test isnan(PiExpTimes(NaN,NaN))
    @test isnan(PiExpTimes(0, Inf))

    @test iszero((0)*Pi)
    @test iszero(PiExpTimes(0, 2))
    @test iszero(PiExpTimes(0, 0))
    @test iszero(Pi^-Inf)
    @test !iszero(PiExpTimes(0, NaN))
end

@testset "show" begin
    @testset "Int" begin
        @test repr(Pi) == "Pi"
        @test repr(2Pi) == "2Pi"
        @test repr(0*Pi) == "0"
    end

    @testset "Rational" begin
        @test repr((1//1)*Pi) == "(1//1)Pi"
        @test repr((2//3)*Pi) == "(2//3)Pi"
        @test repr((2//3)*Pi*Pi) == "(2//3)Pi^2"
    end

    @testset "float" begin
        @test repr(1.0Pi) == "Pi"
        @test repr(2.0Pi) == "2.0*Pi"
        @test repr(0.0Pi) == "0.0"
        @test repr(-0.0Pi) == "-0.0"

        @testset "inf nan" begin
            @test repr(Inf*Pi) == "Inf"
            @test repr(-Inf*Pi) == "-Inf"
            @test repr(NaN*Pi) == "NaN"
        end
    end

    @testset "PiExpTimes" begin
        @test repr(Pi*Pi) == "Pi^2"
        @test repr(Pi*Pi*Pi) == "Pi^3"

        @testset "nested" begin
            @test repr(PiExpTimes{PiExpTimes}(1)) == "Pi^0"
        end
    end

    @testset "complex" begin
        @test repr(im*Pi) == "0 + Pi*im"

        @testset "PiExpTimes" begin
            @test repr(im*Pi*Pi) == "0 + Pi^2*im"
        end

        @testset "float" begin
            @test repr(im*1.2*Pi) == "0.0 + 1.2*Pi*im"
            @test repr((-0.0)*Pi*im) == "-0.0 - 0.0*im"
        end

        @testset "Rational" begin
            @test repr(Pi + im*(2//3)Pi) == "(1//1)Pi + (2//3)Pi*im"
            @test repr(Pi*Pi*(2//3)*im) == "0//1 + (2//3)Pi^2*im"
        end
    end

    @testset "non-finite" begin
        @test repr(Pi^NaN) == repr(NaN)
        @test repr(Pi^Inf) == repr(Inf)
    end
end
