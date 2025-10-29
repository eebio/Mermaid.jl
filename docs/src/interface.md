# Mermaid Interface

In this section, we will:

* See the integrator interface that Mermaid uses,
* Find out how this fits with [CommonSolve](@extref CommonSolve index),

## Interface requirements

The interface used in Mermaid is compatible with, and uses, the SciML [CommonSolve](@extref CommonSolve index) interface.
In particular, each Component is some immutable problem type that stores required data to solve a Component (for example, an [ODEProblem](@extref DiffEq types/ode_types)).
We then also have an Integrator for each Component which stores the current state and handles solving over time.

## Component

There are two categories of Components in Mermaid.
`AbstractTimeDependentComponent`s and `AbstractTimeIndependentComponent`s, both of which are `AbstractComponent`s.

A component must contain the following fields:

* `name` which is a subtype of `AbstractString`,
* `time_step` which is a `Float64`,
* `state_names` which is a `Dict` which maps from variable names in the form of a `String` ([variable name format](@ref "Variable Names")) to some other object (typically a numerical index or symbol) which can be used to identify that variable internally.

## Integrator

A `ComponentIntegrator` is a mutable struct which can be freely modified over a simulation, since it will be typically discarded at the end.
It only stores the current state of that component, and has some associated functions for handling this state.

!!! warning
    An Integrator is required for all Components, even `TimeIndependentComponent`s. They just store the state that they were last called with.

## Functions

Most of the interface for Mermaid is built around functions that manipulate the integrator state.

### step!

The `step!` function should advance the `ComponentIntegrator` one time step (defined by `Component.time_step`).

Its signature should be `step!(comp)`.

### init

The `init` function takes as inputs a `Component` and returns the corresponding `ComponentIntegrator`.

Its function signature should be `init(comp)` and return a subtype of [AbstractComponentIntegrator](@ref).

### getstate and setstate!

The function `getstate` reads the state of an `AbstractComponentIntegrator` at a given `ConnectedVariable` and returns it.
The `setstate!` function similarly mutates the current state of the `ComponentIntegrator` to assign an inputted value to the `ConnectedVariable`.
Additionally, `getstate` and `setstate!` should have methods that do not take a `ConnectedVariable`, and instead return the full state (in a form that should be compatible between the two functions, but need not be documented).

The function signatures should be `getstate(integrator, variable)/setstate!(integrator, variable, value)` and `getstate(integrator)/setstate!(integrator, state)`.

### gettime and settime!

Simply returns the current simulated time of a `ComponentIntegrator`. Defaults to calling `getstate` and `setstate!` with the special variable `"#time"`.

The function signatures should be `gettime(integrator, variable)/settime!(integrator, variable, value)`.

### variables

This function should return all the variables supported by the component, including special variables.

The function signature should be `variables(component)` and return an iterable through strings.

## Variable Names

Mermaid has its own interface for connections too.
A [Connector](@ref) uses instances of [ConnectedVariables](@ref ConnectedVariable) as part of applying the connections.
If you want to write your own components, it is worth learning how the [ConnectedVariables](@ref ConnectedVariable) are defined.

```@docs; canonical=false
ConnectedVariable
```

We can see that each `ConnectedVariable` has two names, the first is the component name and the second is the variable name (which should match `state_names` for that component).

It also has a variable index, which allows you to access only some parts of a full vector. For example, you may index on only some agents in an AgentsComponent.

Finally, it has a duplicated index. Like the variable index, it allows you to act on some subset, but in this case, it is a subset of integrators from a [DuplicatedComponent](@ref).

A `ConnectedVariable` can also be constructed from a single string input, ie `ConnectedVariable("component[duplicatedindex].variable[variableindex]")`.

```@docs; canonical=false
ConnectedVariable(::AbstractString)
```

## Examples

To see some examples, look in the `ext` directory on GitHub.
