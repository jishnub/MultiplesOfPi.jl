# MultiplesOfPi

[![Build Status](https://travis-ci.com/jishnub/MultiplesOfPi.jl.svg?branch=master)](https://travis-ci.com/jishnub/MultiplesOfPi.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jishnub/MultiplesOfPi.jl?svg=true)](https://ci.appveyor.com/project/jishnub/MultiplesOfPi-jl)
[![Codecov](https://codecov.io/gh/jishnub/MultiplesOfPi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/MultiplesOfPi.jl)
[![Coveralls](https://coveralls.io/repos/github/jishnub/MultiplesOfPi.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/MultiplesOfPi.jl?branch=master)

# Introduction

This package exports the type `PiExpTimes{N,T}` that satisfies `PiExpTimes{N,T}(x) = x*π^N`, and the type `PiTimes` that is aliased to `PiExpTimes{1}`. It also provides the constant `Pi` for convenience, defined as `PiTimes(1)`, that behaves like `π` except it produces results with higher accuracy in certain trigonometric and algebraic contexts. In most scenarios the numbers `Pi` and `π` are interchangable.

```julia
julia> Pi^2 == π^2
true
```

It's usually possible, and cleaner, to express mathematical relations in terms of `Pi` instead of the more cumbersome `PiExpTimes`, and is recommended unless it's specifically necessary to do so.

## Rationale

The number `π` is represented as an `Irrational` type in julia, and may be computed to an arbitrary degree of precision. In normal course of events it is converted to a float when it encounters another number, for example `2π` is computed by converting both `2` and `π` to floats and subsequently carrying out a floating-point multiplication. This is lossy, as both `2` and `π` may be represented with arbitrary precision. This package delays the conversion of the `π` to a float, treating it as a common factor in algebraic simplifications. This limits floating-point inaccuracies, especially if the terms multiplying `π` are exactly representable in binary. As an added advantage, it uses `sinpi` and `cospi` wherever possible to avoid having to convert `π` to a float altogether.

# Examples

## Arithmetic

Delaying the conversion of `π` to a float results in satisfying mathematical expressions such as

```julia
julia> (1//3)π + (4//3)π == (5//3)π
false

julia> (1//3)Pi + (4//3)Pi == (5//3)Pi
true
```

We may also simplify algebraic expressions involving powers of `Pi` as 

```julia
julia> (2Pi^2//3) // (4Pi//5)
(5//6)Pi

julia> Pi^-2 / 4Pi^3
0.25*Pi^-5
```

The powers of `Pi` cancel as expected, and `Pi^0` is automatically converted to an ordinary real number wherever possible.

```julia
julia> Pi^2 / Pi^2
1.0
```

Expressions involving `Pi` are automatically promoted to `Complex` as necessary, eg.

```julia
julia> (1+im)Pi^3 / 2Pi^2
0.5*Pi + 0.5*Pi*im

julia> (1+im)Pi^3 / 2Pi^2 * 2/Pi
1.0 + 1.0im
```

## Trigonometric functions

The type `PiTimes` uses `sinpi` and `cospi` under the hood when it is used as an argument to `sin` and `cos`. This results in exact results in several contexts where the inaccuracies arise from floating-point conversions.

```julia
julia> cos(3π/2)
-1.8369701987210297e-16

julia> cos(3Pi/2)
0.0

julia> sin(-π)
-1.2246467991473532e-16

julia> sin(-Pi)
-0.0

julia> tan(π/2)
1.633123935319537e16

julia> tan(Pi/2)
Inf
```

We may compute complex exponential exactly:

```julia
julia> exp(im*π/2)
6.123233995736766e-17 + 1.0im

julia> exp(im*Pi/2)
0.0 + 1.0im

# Euler's identity : exp(iπ) + 1 == 0
julia> exp(im*π) + 1
0.0 + 1.2246467991473532e-16im

julia> exp(im*Pi) + 1
0.0 + 0.0im
```

Hyperbolic functions work as expected:

```julia
# cosh(ix) = cos(x)
# Should be exactly zero for x = π/2
julia> cosh(im*π/2)
6.123233995736766e-17 + 0.0im

julia> cosh(im*Pi/2)
0.0 + 0.0im
```

## Look out

### Type-instability

The type `PiExpTimes{N}` stores the exponent as a type-parameter, therefore exponentiation is not type-stable in general. 

### Floating-point promotion

`PiExpTimes{N}` is promoted to float as soon as it encounters a number other than another `PiExpTimes{N}` in additions and subtractions.

```julia
julia> Pi^2 + Pi
13.011197054679151
```

Arrays of mixed types will involve conversion to the lowest exponent, unless a specific type is specified.

```julia
julia> [Pi,Pi^2]
2-element Array{PiExpTimes{1,Float64},1}:
                   Pi
 3.141592653589793*Pi
 
julia> PiExpTimes[Pi,Pi^2] # this preserves the type of each element, but is not a concrete type
2-element Array{PiExpTimes,1}:
   Pi
 Pi^2
 
julia> PiExpTimes{2}[Pi,Pi^2] # Not a concrete type, but preserves the type of Pi^2 while converting Pi
2-element Array{PiExpTimes{2,T} where T<:Real,1}:
 0.3183098861837907*Pi^2
                    Pi^2
                    
julia> PiExpTimes{2,Float64}[Pi,Pi^2] # Concrete type
2-element Array{PiExpTimes{2,Float64},1}:
 0.3183098861837907*Pi^2
                    Pi^2
```

# Installation

Install the package using 

```julia
pkg> add MultiplesOfPi
```

# Related packages

- [IrrationalExpressions.jl](https://github.com/jishnub/IrrationalExpressions.jl.git)
- [Tau.jl](https://github.com/JuliaMath/Tau.jl)
- [UnitfulAngles.jl](https://github.com/yakir12/UnitfulAngles.jl)
