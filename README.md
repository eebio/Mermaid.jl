# Mermaid.jl

[![Run tests](https://github.com/mjowen/Mermaid.jl/actions/workflows/test.yml/badge.svg)](https://github.com/mjowen/Mermaid.jl/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/mjowen/Mermaid.jl/graph/badge.svg?token=XRLUZB8FQS)](https://codecov.io/gh/mjowen/Mermaid.jl)

Mermaid.jl is a general purpose component-based simulation tool in Julia.
It allows users to connect arbitrary Julia models from a wide range of packages to produce complex coupled simulations.

* Any simulation that can be performed in Julia can be included as a component, from DifferentialEquations.jl, to Agents.jl and beyond.
* Simple interface for extending components to tools not yet supported.
