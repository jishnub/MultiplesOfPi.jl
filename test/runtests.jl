using MultiplesOfPi
using Test

@testset "MultiplesOfPi.jl" begin
    @testset "Constructor" begin
        @test PiTimes(1) ≈ π
        @test PiTimes(1.5) ≈ π*1.5
        @test PiTimes(1//2) ≈ π/2
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
    @testset "algebra" begin
        p1 = PiTimes(3)
        p2 = PiTimes(4)

        @testset "negation" begin
            @test -p1 isa PiTimes{Int}
            @test -p1 == PiTimes(-p1.x)
        end 

        @testset "addition" begin
            @test p1 + p1 == PiTimes(2p1.x)
            @test p1 + p2 == PiTimes(p1.x + p2.x)
        end

        @testset "subtraction" begin
            @test p1 - p1 == PiTimes(0)
            @test p1 - p2 == PiTimes(p1.x - p2.x)
        end

        @testset "multiplication" begin
            @test p1*2 == PiTimes(p1.x*2)
            @test 2*p1 == PiTimes(p1.x*2)
            @test 2p1 == PiTimes(p1.x*2)
            @test p1*p1 ≈ (p1.x*π)^2
            @test p1*p2 ≈ (p1.x*p2.x)*π^2
            @test p1 * im == PiTimes(0) + PiTimes(p1.x)*im 
            @test im* p1 == PiTimes(0) + PiTimes(p1.x)*im
            @test (1+im)*p1 == PiTimes(p1.x) + PiTimes(p1.x)*im
        end

        @testset "division" begin
            @test p1/p1 == 1
            @test p1/p2 == p1.x/p2.x == 3/4 # works as 1/4 can be stored exactly
            @test p2/p1 == p2.x/p1.x

            @test p1/π === p1.x
            @test PiTimes(7.3)/π === 7.3
            @test π/PiTimes(1/7.3) ≈ 7.3 # within rounding errors
            @test π/PiTimes(1//7) === 7//1
        end
    end
    @testset "trigonometric functions" begin
        x = 1.5
        @test sin(PiTimes(x)) == sinpi(x) ≈ sin(π*x)
        x = 1
        @test cos(PiTimes(x)) == cospi(x) ≈ cos(π*x)
        x = 1.1
        @test cis(PiTimes(x)) == cospi(x) + im*sinpi(x) ≈ cos(π*x) + im*sin(π*x)

        a = sincos(PiTimes(x))
        @test a[1] == sinpi(x)
        @test a[2] == cospi(x)

        @testset "exp" begin
            @test exp(im*PiTimes(0)) == cis(PiTimes(0)) == 1
            @test exp(im*PiTimes(1/2)) == cis(PiTimes(1/2)) == im
            @test exp(im*PiTimes(1)) == cis(PiTimes(1))== -1
            @test exp(im*PiTimes(3/2)) == cis(PiTimes(3/2)) == -im
        end
    end
    @testset "hyperbolic" begin
        @testset "cosh" begin
            @test cosh(im*PiTimes(1/2)) == 0
        end
        @testset "sinh" begin
            @test sinh(im*PiTimes(-1)) == 0
        end
    end
end
