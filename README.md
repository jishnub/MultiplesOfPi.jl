# MultiplesOfPi

[![Build Status](https://travis-ci.com/jishnub/MultiplesOfPi.jl.svg?branch=master)](https://travis-ci.com/jishnub/MultiplesOfPi.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jishnub/MultiplesOfPi.jl?svg=true)](https://ci.appveyor.com/project/jishnub/MultiplesOfPi-jl)
[![Codecov](https://codecov.io/gh/jishnub/MultiplesOfPi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/MultiplesOfPi.jl)
[![Coveralls](https://coveralls.io/repos/github/jishnub/MultiplesOfPi.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/MultiplesOfPi.jl?branch=master)

# Introduction

This package introduces the type `PiTimes` that automatically uses the functions `sinpi` and `cospi` instead of `sin` and `cos` to produce accurate results by avoiding floating-point inaccuracies to some extent. It also provides the constant `Pi` for convenience, defined as `PiTimes(1)`, which behaves like `pi` except it produces results with higher accuracy in certain trigonometric contexts.

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

## Trigonometric functions

The `PiTimes` type uses `sinpi` and `cospi` under the hood when it is used as an argument to `sin` and `cos`. This results in exact results in several contexts where the inaccuracies arise from floating-point conversions.

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

It automatically promotes to `Complex` as necessary, so we may compute complex exponentials exactly for some arguments:

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

# Installation

Install the package using 

```julia
pkg> add MultiplesOfPi
```

# Related packages

- [IrrationalExpressions.jl](https://github.com/jishnub/IrrationalExpressions.jl.git)