# Mermaid.jl

![Mermaid.jl](https://raw.githubusercontent.com/eebio/Mermaid.jl/refs/heads/registry-prep/docs/src/assets/logo-full.svg)
[![Run tests](https://github.com/mjowen/Mermaid.jl/actions/workflows/test.yml/badge.svg)](https://github.com/mjowen/Mermaid.jl/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/mjowen/Mermaid.jl/graph/badge.svg?token=XRLUZB8FQS)](https://codecov.io/gh/mjowen/Mermaid.jl)

Mermaid.jl is a hybrid and multiscale simulation enviroment in Julia.

Complex simulations can be produced by connecting together components from a wide range of Julia modeling tools.

Its key features are:

1. It is particularly well suited towards hybrid (continuous and discrete time), multiscale and nested systems. With direct support for nesting models within other models (Agent-based models where each agent solves an ODE, for example).
2. Models can be specified as arbitrary Julia code, include calls to other programming languages such as C or Python.
3. Out-of-the-box support for Agents.jl, DifferentialEquations.jl (and related packages), Surrogates.jl, and MethodOfLines.jl.
4. A simple integrator interface allows extensions to previously unsupported components.
