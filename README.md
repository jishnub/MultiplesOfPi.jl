# MultiplesOfPi

[![Build Status](https://travis-ci.com/jishnub/MultiplesOfPi.jl.svg?branch=master)](https://travis-ci.com/jishnub/MultiplesOfPi.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jishnub/MultiplesOfPi.jl?svg=true)](https://ci.appveyor.com/project/jishnub/MultiplesOfPi-jl)
[![Codecov](https://codecov.io/gh/jishnub/MultiplesOfPi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/MultiplesOfPi.jl)
[![Coveralls](https://coveralls.io/repos/github/jishnub/MultiplesOfPi.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/MultiplesOfPi.jl?branch=master)

# Introduction

This package introduces the type `PiTimes` that automatically uses the functions `sinpi` and `cospi` instead of `sin` and `cos` to produce accurate results by avoiding floating-point inaccuracies to some extent. It also provides the constant `Pi` for convenience, defined as `PiTimes(1)`, which behaves like `pi` except it produces results with higher accuracy in certain trigonometric contexts.

For example:

```julia
julia> cos(3π/2)
-1.8369701987210297e-16

julia> cos(3Pi/2)
0.0

julia> sin(-π)
-1.2246467991473532e-16

julia> sin(-Pi)
-0.0
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

The function `exp(im*x)` for a finite real `x` might also be computed using `cis`, as
```julia
julia> cis(Pi/2)
0.0 + 1.0im
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
pkg> add https://github.com/jishnub/MultiplesOfPi.jl.git

julia> using MultiplesOfPi
```
# To-do

Add support for `tan`