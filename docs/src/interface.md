# Mermaid Interface

In this section, we will:

* See the integrator interface that Mermaid uses,
* Find out how this fits with [CommonSolve](@extref CommonSolve index),
* Look at how the [AgentsComponent](@ref) is defined, as an example.

## Interface requirements

The interface used is Mermaid is compatible with, and uses, the SciML [CommonSolve](@extref CommonSolve index) interface.
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
It typically only stores the current state of that component, and has some associated functions for handling this state.

In addition to having some functions defined on the integrator, we also have some required fields.

* `outputs` and `inputs` are `Dict`s for mapping from a [ConnectedVariable](@ref) to its current value. The values in these `Dict`s are set by Mermaid but the initial creation of this `Dict` (in particular, its keys) need to be defined in the `init` function. The values from `inputs` should be used to set the state during the `step!` function.
* The original `Component` should also be stored under the field name `component`. This allows access to any problem settings, or other fields like the `name`.

!!! todo "setstate! and inputs"
    Couldn't the state be set by mermaid, rather than in the step function?

!!! warning
    An Integrator is required for all Components, even `TimeIndependentComponent`s. They just store the state that they were last called with.

## Functions

### step!

The `step!` function should advance the `ComponentIntegrator` (its only input) one time step (defined by `Component.time_step`).
This should also include reading the `inputs` field and updating the integrators internal state accordingly.

### init

The `init` function takes as inputs a `Component` and a `Vector{Connector}` and returns the corresponding `ComponentIntegrator`.
This includes the process of creating the `inputs` and `outputs` `Dict`s which require the keys to be specified, i.e. the `init` function should read through the `Vector{Connector}` and look for any outputs or inputs into this component, and add the `ConnectedVariable`s as keys to the `inputs` and `outputs` `Dict`s.

### getstate and setstate!

The function `getstate` function reads the state of a `ComponentIntegrator` at a given `ConnectedVariable` and returns it.
The `setstate!` function similarly mutates the current state of the `ComponentIntegrator` to assign an inputted value to the `ConnectedVariable`.

### gettime

Simply returns the current simulated time of a `ComponentIntegrator`.

## AgentsComponent Example

!!! todo "Code example"
    We should include the code of the Agents component here too. probably with literate

```@docs; canonical=false
AgentsComponent
Mermaid.AgentsComponentIntegrator
Mermaid.step!(::Mermaid.AgentsComponentIntegrator)
CommonSolve.init(::Mermaid.AgentsComponent, ::Vector{Connector})
Mermaid.getstate(::Mermaid.AgentsComponentIntegrator, ::Mermaid.ConnectedVariable)
Mermaid.setstate!(::Mermaid.AgentsComponentIntegrator, ::Mermaid.ConnectedVariable, value)
Mermaid.gettime(::Mermaid.AgentsComponentIntegrator)
```

## Variable Names
