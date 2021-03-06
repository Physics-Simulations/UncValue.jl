![tests](https://github.com/Physics-Simulations/UncValue.jl/workflows/tests/badge.svg) [![codecov](https://codecov.io/gh/Physics-Simulations/UncValue.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Physics-Simulations/UncValue.jl) [![version](https://juliahub.com/docs/UncValue/version.svg)](https://juliahub.com/ui/Packages/UncValue/2sLzV)

# Uncertainty-Value
Simple class to evaluate the uncertainty for complex or very long calculations given the initial values together with its uncertainty.

# How-To

## Installation
The recommended way to install UncValue in your computer is via `Pkg` by writting
```julia
using Pkg
Pkg.add("UncValue")
```

Alternatively, you can download the [release files](https://github.com/Physics-Simulations/UncValue.jl/releases) and install it manually.


## Usage
The way it works is simple, first import the script as
```julia
using UncValue
```
then initialise your `Value` variables (numbers, lists, matrices...) as
```julia
pi = Value(3.14159, 0.00011) # number variable 3.14159 +/- 0.00011
A = [pi; Value(2.718, 0.036); Value(1.61803398875, 29e-11)] # array with 3 elements
M = Value(rand(3,5), rand(3,5)*0.056) # 3x5 matrix
```

- `pi` is just a number variable with uncertainty
- `A` is a list of values, each one with each own uncertainty
- `M` is a 3x5 value matrix (not a matrix of values) where the uncertainty is separated from the value, so this class only works as a container for keeping them together but some operations will not work properly (like multiplication). To initialize the matrix of values correctly we should do it as the list.

Perform any operation you want between Value(s):
- Binary operators: `+`, `-`, `*`, `/`, `\`, `^`...
- Unary operators: `abs`, `exp` (base 2, 10 and general), `log` (base 2, 10 and general), `sqrt`, `cbrt`, trigonometric and inverse function, hyperbolic and inverse functions...
- Comparison: `>=`, `>`, `==`, `!=`, `<`, `<=`...

A complete list of compateble operations can be found in the [Julia documentation](https://docs.julialang.org/en/v1/base/math/).

# Contributors
- [Adrià Labay Mora](https://labay11.github.io/)
- [Àlex Giménez Romero](https://github.com/agimenezromero)

# License
      Copyright 2020 Physics-Simulations

      Licensed under the Apache License, Version 2.0 (the "License");
      you may not use this file except in compliance with the License.
      You may obtain a copy of the License at

          http://www.apache.org/licenses/LICENSE-2.0

      Unless required by applicable law or agreed to in writing, software
      distributed under the License is distributed on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
      See the License for the specific language governing permissions and
      limitations under the License.
