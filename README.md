# MultiplesOfPi

[![Build Status](https://travis-ci.com/jishnub/MultiplesOfPi.jl.svg?branch=master)](https://travis-ci.com/jishnub/MultiplesOfPi.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jishnub/MultiplesOfPi.jl?svg=true)](https://ci.appveyor.com/project/jishnub/MultiplesOfPi-jl)
[![Codecov](https://codecov.io/gh/jishnub/MultiplesOfPi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jishnub/MultiplesOfPi.jl)
[![Coveralls](https://coveralls.io/repos/github/jishnub/MultiplesOfPi.jl/badge.svg?branch=master)](https://coveralls.io/github/jishnub/MultiplesOfPi.jl?branch=master)

# Introduction

This package introduces the type `PiTimes` that automatically uses the functions `sinpi` and `cospi` instead of `sin` and `cos` to produce accurate results by avoiding floating-point inaccuracies to some extent. For example:

```julia
julia> cos(3π/2)
-1.8369701987210297e-16

julia> cos(PiTimes(3/2))
0.0

julia> sin(-π)
-1.2246467991473532e-16

julia> sin(PiTimes(-1))
-0.0
```

# Installation

Install the package using 

```julia
pkg> add https://github.com/jishnub/MultiplesOfPi.jl.git

julia> using MultiplesOfPi
```

# To-do

Add support for complex numbers
