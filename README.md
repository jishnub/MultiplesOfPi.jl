# MultiplesOfPi

[![CI](https://github.com/jishnub/MultiplesOfPi.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/jishnub/MultiplesOfPi.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/jishnub/MultiplesOfPi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/MultiplesOfPi.jl)
[![Coverage Status](https://coveralls.io/repos/github/jishnub/MultiplesOfPi.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/MultiplesOfPi.jl?branch=master)
[![docs:stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jishnub.github.io/MultiplesOfPi.jl/stable)
[![docs:dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jishnub.github.io/MultiplesOfPi.jl/dev)

# Introduction

This package exports the type `PiExpTimes` that satisfies `PiExpTimes(x, n) = x*π^n`. It also provides the constant `Pi` for convenience, defined as `PiExpTimes(1, 1)`, that behaves like `π` except it produces results with higher accuracy in certain trigonometric and algebraic contexts. In most scenarios the numbers `Pi` and `π` are interchangable. Expressing mathematical relations in terms of `Pi` instead of the more cumbersome `PiExpTimes` is usually cleaner, and is recommended unless it's specifically necessary to do otherwise.

## Rationale

The number `π` is represented as an `Irrational` type in julia, and may be computed to an arbitrary degree of precision. In normal course of events it is converted to a float when it encounters another number, for example `2π` is computed by converting both `2` and `π` to floats and subsequently carrying out a floating-point multiplication. This is lossy, as both `2` and `π` may be represented with arbitrary precision. This package delays the conversion of the `π` to a float, treating it as a common factor in algebraic simplifications. This limits floating-point inaccuracies, especially if the terms multiplying `π` are exactly representable in binary. As an added advantage, it uses `sinpi` and `cospi` wherever possible to avoid having to convert `π` to a float altogether.

# Features

## Arithmetic

Delaying the conversion of `π` to a float results in precise mathematical expressions such as

```julia
julia> (1//3)π + (4//3)π == (5//3)π
false

julia> (1//3)Pi + (4//3)Pi == (5//3)Pi
true

julia> float(sqrt(pi)^2) == float(pi)
false

julia> float(sqrt(Pi)^2) == float(Pi)
true
```

We may also simplify algebraic expressions involving powers of `Pi` as

```julia
julia> (2Pi^2//3) // (4Pi//5)
(5//6)Pi

julia> Pi^-2 / 4Pi^3
0.25*Pi^-5
```

Expressions involving `Pi` are automatically promoted to `Complex` as necessary, eg.

```julia
julia> (1+im)Pi^3 / 2Pi^2
0.5*Pi + 0.5*Pi*im

julia> (1+im)Pi^3 / 2Pi^2 * 2/Pi
Pi^0 + Pi^0*im
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

### Interactions with π

The irrational number `π` is usually aggressively converted to `Pi`, eg:

```julia
julia> π * Pi
Pi^2
```

This ensures that subsequent calculation would not get promoted to a floating-point type. However if this behavior is not desired then one may specify the type explicitly while constructing the object as

```julia
julia> PiExpTimes{Irrational{:π}}(π)
π*Pi^0
```

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

Note that addition and subtraction are not type-stable at present by design.

# Performance

In general using `Pi` instead of `pi` will be less performant, as `pi` is aggressively promoted to a floating-point value in most calculations. The use of `Pi` is mainly intended for tests where exact fractions of `Pi` are desirable.

# Installation

Install the package using

```julia
pkg> add MultiplesOfPi
```

# Related packages

- [IrrationalExpressions.jl](https://github.com/jishnub/IrrationalExpressions.jl.git)
- [Tau.jl](https://github.com/JuliaMath/Tau.jl)
- [UnitfulAngles.jl](https://github.com/yakir12/UnitfulAngles.jl)
