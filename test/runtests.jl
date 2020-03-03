using MultiplesOfPi
using Test

@testset "MultiplesOfPi.jl" begin
    @testset "Constructor" begin
        @test PiTimes(1) ≈ π
        @test PiTimes(1.5) ≈ π*1.5
        @test PiTimes(1//2) ≈ π/2
    end
    @testset "Ordering" begin
        @test PiTimes(1) < 4
        @test PiTimes(2) > 4
        @test PiTimes(2) > π
        @test PiTimes(2) < PiTimes(4)
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
    end
end
