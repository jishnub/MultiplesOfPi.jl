# MultiplesOfPi

[![Build Status](https://travis-ci.com/jishnub/MultiplesOfPi.jl.svg?branch=master)](https://travis-ci.com/jishnub/MultiplesOfPi.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jishnub/MultiplesOfPi.jl?svg=true)](https://ci.appveyor.com/project/jishnub/MultiplesOfPi-jl)
[![Codecov](https://codecov.io/gh/jishnub/MultiplesOfPi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/MultiplesOfPi.jl)
[![Coverage Status](https://coveralls.io/repos/github/jishnub/MultiplesOfPi.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/MultiplesOfPi.jl?branch=master)

# Introduction

This package exports the type `PiExpTimes{N,T}` that satisfies `PiExpTimes{N,T}(x) = x*π^N`, and the type `PiTimes` that is aliased to `PiExpTimes{1}`. It also provides the constant `Pi` for convenience, defined as `PiTimes(1)`, that behaves like `π` except it produces results with higher accuracy in certain trigonometric and algebraic contexts. In most scenarios the numbers `Pi` and `π` are interchangable.

```julia
julia> Pi^2 == π^2
true
```

Expressing mathematical relations in terms of `Pi` instead of the more cumbersome `PiExpTimes` is usually cleaner, and is recommended unless it's specifically necessary to do otherwise. One such scenario is exponentiation, more on which is presented below.

## Rationale

The number `π` is represented as an `Irrational` type in julia, and may be computed to an arbitrary degree of precision. In normal course of events it is converted to a float when it encounters another number, for example `2π` is computed by converting both `2` and `π` to floats and subsequently carrying out a floating-point multiplication. This is lossy, as both `2` and `π` may be represented with arbitrary precision. This package delays the conversion of the `π` to a float, treating it as a common factor in algebraic simplifications. This limits floating-point inaccuracies, especially if the terms multiplying `π` are exactly representable in binary. As an added advantage, it uses `sinpi` and `cospi` wherever possible to avoid having to convert `π` to a float altogether.

# Features

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

## Algebra

We may convert between types having different exponents without losing accuracy.

```julia
julia> convert(PiExpTimes{3},Pi)
Pi^-2*Pi^3

julia> convert(PiExpTimes{3},Pi) == Pi
true
```

Such an expression may be reduced to a simple form using `simplify`.

```julia
julia> convert(PiExpTimes{3},Pi) |> MultiplesOfPi.simplify
Pi
```

## Look out

### Type-instability

The type `PiExpTimes{N}` stores the exponent as a type-parameter, therefore exponentiation is not type-stable in general.

### Constructor-abuse to avoid nesting

`PiExpTimes{N}(PiExpTimes{M})` is automatically simplified to `PiExpTimes{N+M}`. This is an abuse of Julia's constructors as the type of the object changes, however this avoids nested expressions that have performance issues.

```julia
julia> PiExpTimes{2}(PiExpTimes{3}(4))
4Pi^5

julia> PiExpTimes{2}(PiExpTimes{3}(4)) |> typeof
PiExpTimes{5,Int64}
```

### Interactions with π

The irrational number `π` is usually aggressively converted to `PiTimes(1)`, eg:

```julia
julia> PiTimes(π)
Pi^2
```

This ensures that subsequent calculation would not get promoted to a floating-point type. However if this behavior is not desired then one may specify the type explicitly while constructing the object as 

```julia
julia> PiTimes{Irrational{:π}}(π)
π = 3.1415926535897...*Pi
```

However, it is not possible to convert a number to the type `PiTimes{Irrational{:π}}`.

### Floating-point promotion

Addition and subtraction involving mixed exponents of `Pi` will involve floating-point conversions, and the resulting expression will have the minimum exponent out of the terms being summed.

```julia
julia> Pi + 3Pi
4Pi

julia> Pi + 3Pi^2
10.42477796076938*Pi
```

This fits with the intuition of the expression being factorized as `Pi + 3Pi^2 == Pi*(1 + 3Pi)`.

Note that `π` is promoted to `Pi` in such operations, so we obtain 

```julia
julia> Pi + π
2Pi
```

### Conversion vs construction

Constructors for `PiExpTimes` act as a wrapper and not as a conversion. Conversion to the type retains the value of a number, whereas construction implies multiplication by an exponent of `Pi`.

```julia
julia> PiTimes(1)
Pi

julia> convert(PiTimes,1)
Pi^-1*Pi
```

### Promotion of mixed types

`promote` and `promote_type` work differently with types having different exponents. `promote` converts both the types to one that has the minimum exponent, whereas `promote_type` leaves the exponent as a free parameter.

```julia
julia> promote(Pi,Pi^2) |> typeof
Tuple{PiExpTimes{1,Real},PiExpTimes{1,Real}}

julia> promote_type(typeof(Pi),typeof(Pi^2))
PiExpTimes{N,Int64} where N
```

This is so that structs of `PiTimes` &mdash; such as complex numbers &mdash; do not lose accuracy in conversion. The behaviour is explained below with some examples:

#### Arrays of mixed types

Storing different exponents of `Pi` in an array will in general lead to conversion to a supertype that can store all the values.

```julia
julia> [Pi,Pi^2] # element type is not concrete
2-element Array{PiExpTimes{N,Int64} where N,1}:
   Pi
 Pi^2
```

Such an array will not be the most performant, as the element type is not concrete. This may be avoided by specifying a type while creating the array. A concrete type will not be able to store all the numbers losslessly.

```julia
julia> PiExpTimes{2,Real}[Pi,Pi^2] # exponent is fixed but Real is not a concrete type
2-element Array{PiExpTimes{2,Real},1}:
 Pi^-1*Pi^2
       Pi^2

julia> PiExpTimes{2,Float64}[Pi,Pi^2] # eltype is concrete, but result loses accuracy
2-element Array{PiExpTimes{2,Float64},1}:
 0.3183098861837907*Pi^2
                    Pi^2
```

#### Complex numbers

Complex numbers rely on `promote` to generate the element type

```julia
julia> Complex(Pi,Pi^2)
Pi + Pi*Pi*im

julia> Complex(Pi,Pi^2) |> typeof
Complex{PiExpTimes{1,Real}}
```

In this case converting either the real or imaginary part to a floating-point type would have resulted in a loss of accuracy. Such a type might not be performant, so if a conversion is desired then it might be enforced by specifying the element type while constructing the `Complex` struct:

```julia
julia> Complex{PiTimes{Float64}}(Pi,Pi^2)
Pi + 3.141592653589793*Pi*im

julia> Complex{PiTimes{Float64}}(Pi,Pi^2) |> typeof
Complex{PiExpTimes{1,Float64}}
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
