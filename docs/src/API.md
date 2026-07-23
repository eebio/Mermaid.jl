# API

## Contents

```@contents
Pages = ["API.md"]
Depth = 2:3
```

This page documents the public API and the execution rules used by Mermaid.
The [Tutorial](@ref) is a narrative introduction; this page is a reference for constructing problems, writing connections, implementing components, and reading solutions.

## Core concepts

A Mermaid simulation has four layers:

1. A **component** is an immutable description of a model and its configuration.
2. A **component integrator** is the mutable runtime state created from a component.
3. A [`MermaidProblem`](@ref) stores components, connectors, the global time span, and component time scales.
4. A [`MermaidSolution`](@ref) stores selected variables at selected global times.

`init(component)` creates a component integrator.
`init(problem, algorithm)` creates a [`MermaidIntegrator`](@ref), and `solve!` advances it until the end of the problem's time span.
The convenience function `solve(problem, algorithm; kwargs...)` is provided by CommonSolve and performs initialization followed by solving.

One of the main goals of Mermaid is to make it simple and accessible to expand the functionality to user-defined components.
This is done through the component interface, a set of functions (outlined below) that can be implemented to support all Mermaid functionality and include the component in hybrid simulations.

```julia
init(component)
step!(component_integrator)
getstate(component_integrator)
getstate(component_integrator, variable)
setstate!(component_integrator, state)
setstate!(component_integrator, variable, value)
variables(component_or_integrator)
name(component_or_integrator)
timestep(component_or_integrator)
```

## Component Interface

`name(component)` returns the unique component name.
`timestep(component)` returns the proposed local step size.
`variables(component)` returns the names accepted by `getstate` and `setstate!`, including supported special variables.
Methods on an integrator normally delegate to its component.

```@docs
name
timestep
variables
getstate
setstate!
gettime
settime!
step!(::AbstractComponentIntegrator)
init(::AbstractComponent)
```

`getstate(integrator)` and `setstate!(integrator, state)` operate on the complete state (except for time), the specifics of the return value may optionally be documented.
If the return values are undocumented, then the only requirement is that `setstate!(getstate(int))` doesn't alter the state of the integrator.
The keyed forms operate on one named or special variable.
`getstate(...; copy=true)` returns a deep copy, which is important when retaining mutable state such as an Agents.jl model.

While there are quite a few functions here that need to be defined, the majority of them have default behaviour that you can opt in to.

* `name(integrator)`, `timestep(integrator)`, and `variables(integrator)` can use the default behavior of `name(integrator) = name(integrator.component)` (and other equivalents) if the component is stored as a property or field of the integrator.
* `timestep(component)` and `name(component)` can use the default behaviour of `component.timestep` and `component.name`, if the property or field exists.
* `gettime` and `settime!` will default to attempting to access the `#time` special variable of `getstate` and `setstate!`.
* `getstate` and `setstate!` do not need to implement any of the keyword argument features.

## Problems and algorithms

### `MermaidProblem`

```@docs
MermaidProblem
```

`components` is ordered and component names must be unique because connectors refer to components by name.
The order of `components` determines the order in which the components are stepped when multiple components can be stepped together, although unintential side effects of stepping components is unlikely.
`connectors` is also ordered; that order is significant when one connector updates a value used by a later connector.
The order of `connectors` is much more relevant.
While in most situations, applying a connection does not perturb other variables, there are some exceptions to this.
If special variables, like `#state` or `#integrator` are perturbed, it is possible for the order of connections to be influential.
Similarly, setting `#ids` in a duplicated component will immediately change the size (and possibly order) of its state vector `#states`.
`tspan` is the global time span, marking the start time and end time of the simulation.
Timescales, discussed in [Timestepping](@ref), allow you to specify models using different time units, and relate them together into a single global simulation timescale.

### `MermaidIntegrator`

```@docs
MermaidIntegrator
```

The integrator contains the runtime component integrators, connector list, algorithm, saving configuration, and global current time.
It is mutable and is normally consumed by `solve!`; use `step!(integrator)` for manual advancement.

### Initialization and solving

```@docs
init(::AbstractMermaidProblem, ::AbstractMermaidSolver)
solve!
step!(::AbstractMermaidIntegrator)
```

`solve!` records the initial state when it satisfies `saveat`, then repeatedly calls `step!` until the global end time.
`step!(mermaid_integrator)` dispatches to the algorithm stored in `mermaid_integrator.alg`.

The currently provided algorithm is:

```@docs
MinimumTimeStepper
```

`MinimumTimeStepper` advances the global clock to the earliest next component event, subject to a possible earlier requested `saveat` time.
It applies eligible connectors before stepping components, then steps every component whose next local event reaches the new global time.
The specifics of this algorithm are documented in [Timestepping](@ref).

## Connections

### Connected variables

```@docs
ConnectedVariable
ConnectedVariable(::AbstractString)
fullname
```

A connected variable has four pieces:

| Syntax/Fullname | Meaning |
| --- | --- |
| `component.variable` | A variable in a component |
| `component.variable[i]` | An index into a variable |
| `component[j].variable` | Instance `j` of a duplicated component |
| `component[j].variable[i]` | Both kinds of indexing |

The component and variable names are strings.
A variable index is an integer or an integer collection/range.
A duplicated index selects IDs in a [`DuplicatedComponent`](@ref), or rather there corresponding `#state`.
String construction parses the index expressions as Julia expressions, so use literal integer indices and ranges such as `1`, `1:5`, or `[1, 3]`.

`fullname(variable)` produces the canonical string representation.

### Connectors

```@docs
Connector
runconnection
runconnection!
```

The usual constructor is:

```julia
Connector(; inputs=["source.x"], outputs=["target.y"], func=nothing)
```

Connector semantics are:

* Inputs are read in the order given by `inputs`.
* With `func=nothing`, one input is passed unchanged to every output. With multiple inputs, the values are collected into a vector and that vector is passed to every output.
* With a function, Mermaid calls `func(input1, input2, ...)`; the function receives positional arguments, not the input vector. Its return value is passed unchanged to every output.
* `runconnection` computes and returns the output without mutating state. `runconnection!` computes the output and calls `setstate!` for each output.
* Outputs are written in the order given by `outputs`. A connector does not perform automatic unit conversion, interpolation, broadcasting, reshaping, or spatial mapping. The function must return a value compatible with every output. Mapping between spaces or shapes must be explicit in `func`.
* A connector with no outputs is valid and is useful for side effects such as recording an animation frame. Its function's return value is ignored.

Connectors are applied in `MermaidProblem.connectors` order.
During a `MinimumTimeStepper` event, a connector is applied only when every input is no later than every output in global time.
Eligible connectors are applied before component steps.
Later connectors can therefore observe changes made by earlier connectors, but a connector does not automatically run again after the component step unless a later event makes it eligible.

### Connection mapping functions

Variable indices are passed to the component's `getstate`/`setstate!` implementation.
Their exact shape follows the integration.
For a duplicated component, the duplicated index selects instances and the variable index selects a part of each instance's state.
Both indices may be used together.

Mermaid has no separate built-in spatial-map abstraction.
A mapping between coordinate systems, grids, resolutions, units, or state layouts must be implemented as the connector function:

```julia
map_grid_to_agents(grid_values, positions) = interpolate(grid_values, positions)
Connector(inputs=["pde.concentration", "agents.positions"],
          outputs=["agents.concentration"], func=map_grid_to_agents)
```

## Timestepping

`MermaidProblem` accepts one `timescales` value per component:

```julia
MermaidProblem(components=components, connectors=connectors,
               tspan=(0.0, 100.0), timescales=[1.0, 0.1])
```

For component `i`, Mermaid time is:

$$t_{global,i} = s_i t_{local,i},$$

where $s_i$ is `timescales[i]`.
A component with `s_i = 0.1` advances ten local time units for every one global time unit.

The minimum stepper chooses the smallest upcoming value of `(local_time + timestep) * timescale`.
It applies connectors whose input global times are no later than their output global times and steps all components whose next event reaches that global time.
Some minor changes to the time values may occur to minigate floating-point roundoff accumulation.

This method has a nice property compared with other possible methods.
Notably, it guarantees that possible connection applications cannot be jumped over.
If the method synchronised at different timepoints, it is possible that a connection that could have been applied if time was treated continuously, would be ignored.

Connections use the most recently available state at the synchronization event.
They do not interpolate between component states and do not guarantee identical local times.

## Saving and solutions

### Saving options

Saving is configured during initialization:

```julia
init(problem, alg; save_vars=nothing, saveat=nothing)
```

`save_vars` accepts:

* `nothing` (default): all non-special variables returned by `variables`.
* `:all`: all variables, including special variables beginning with `#`.
* `:none` or `String[]`: no variables, while save times are still recorded.
* A vector of strings: exactly the listed connected-variable names. Indices may select subsets such as `"forest.life[1]"` or `"tree[1:10].life"`.

`saveat` accepts:

* `nothing` (default): save after initialization and after every Mermaid event.
* A number `Δt`: save at `tspan[1]:Δt:tspan[2]`.
* A vector of times to stop and save at.
* A function `(integrator, t) -> Bool`: save whenever it returns `true`. This is simply checked at each of the timepoints the algorithm was already going to stop at, and checks this function. No attempt is made to ensure these `saveat` points are not missed.

### `MermaidSolution`

```@docs
MermaidSolution
```

A solution has `sol.t`, the saved global times, and `sol.u`, a dictionary whose keys are `ConnectedVariable`s.
Each dictionary value has one entry per saved time.

```julia
sol["component.variable"]
sol[ConnectedVariable("component.variable")]
sol[3]
sol(2.5)
```

String indexing is converted to a `ConnectedVariable` and returns the saved history.
`ConnectedVariable` indexing is equivalent and is useful for programmatic keys.
Integer indexing returns a new [`MermaidSolution`](@ref) containing only the data at the indexed position for each variable.

Calling `sol(t)` returns a one-time-point solution from linear interpolation between surrounding saved states.
The time must be within `sol.t`.
Interpolation applies component-wise arithmetic to saved values, so it is appropriate for numeric states but generally not meaningful for non-numeric objects such as an Agents.jl model.

!!! todo
    Numeric data should be linearly interpolated, but anything else should use a constant interpolation

## Components

All time-dependent components expose `name`, `timestep`, `variables`, `init`, `step!`, `getstate`, `setstate!`, `gettime`, and `settime!` through the common interface.
The following sections document component-specific construction and special variables.

### Differential Equations Components

```@docs
DEComponent
DEComponentIntegrator
```

`DEComponent(model, alg; name, timestep, intkwargs, state_names)` wraps a SciML DifferentialEquations problem.
`state_names` maps strings to state-vector indices or symbolic variables.

The special variables are:

* `#time`, `integrator.t`;
* `#state`, `integrator.u`;
* `#integrator`, `integrator`;

`getstate` returns the `#state`.

### Agents Components

```@docs
AgentsComponent
AgentsComponentIntegrator
```

`AgentsComponent(model; name, state_names, timestep)` wraps an Agents.jl `StandardABM`.
A `state_names` value may refer to a model property or an agent property. Without a variable index, an agent property returns values for all agents; with an index, it selects agent IDs.
If the property is indexable, please use the `func` field of the connector to access specific property indexes, rather than the `variableindex` field which is reserved for agent IDs.

The special variables are:

* `#time`, a separate clock to `abmtime(model)`;
* `#model`, the current `StandardABM`;
* `#ids`, the vector of all current agent IDs (cannot be used with `setstate!`).

`getstate` returns the `#model`.

### MethodOfLines Components

```@docs
MOLComponent
MOLComponentIntegrator
```

`MOLComponent` wraps a MethodOfLines-generated SciML problem. `state_names` identifies entries or groups of entries in the discretized state.

The special varialbes are:

* `#time`, `integrator.t`;
* `#state`, `integrator.u`;
* `#integrator`, `integrator`;

`getstate` returns the `#state`.

### JumpProcesses Components

```@docs
JumpComponent
JumpComponentIntegrator
```

`JumpComponent` wraps a JumpProcesses `JumpProblem`.

The special variables are:

* `#time`, `integrator.t`;
* `#state`, `integrator.u`;
* `#integrator`, `integrator`;

`getstate` returns the `#state`.

### TrixiParticles Components

```@docs
TrixiParticlesComponent
TrixiParticlesComponentIntegrator
```

`TrixiParticlesComponent(semi, alg; name, timestep, intkwargs, tspan, state_names)` semidiscretizes the supplied TrixiParticles object and wraps the resulting problem.

The special variables are:

* `#time`, `integrator.t`;
* `#state`, `integrator.u`;
* `#integrator`, `integrator`;
* `#semi`, the semidiscretization (cannot be used with `setstate!`).

`getstate` returns the `#state`.

### Surrogate Components

```@docs
SurrogateComponent
SurrogateComponentIntegrator
```

`SurrogateComponent(component, surrogate, lower_bound, upper_bound; kwargs...)` trains a surrogate to approximate one component `step!`.
Training occurs during initialization using Sobol samples within the supplied bounds.
The surrogate currently does not add explicit time as an input, so non-autonomous source components are not recommended.

### Duplicated Components

```@docs
DuplicatedComponent
DuplicatedComponentIntegrator
```

`DuplicatedComponent(component, init_states; instances, name, timestep, state_names, default_state)` reuses one component definition for many independent states.
With a fixed `instances` value, positions initialize as `1:instances`.
With `instances=nothing`, the component is flexible and initially has no IDs; connect its `#ids` variable to create and remove states.

The special variables are:

* `#time`;
* `#ids`, current integer IDs in current state order;
* `#states`, complete vector of instance states;
* `#init_states`, initial states used when an ID is first created.

Any special variables of the underlying system (other than the ones defined above) are also available.

When `#ids` is set, it is compared against the previous value of `#ids`, existing IDs retain their states, missing IDs are removed, and new IDs receive a copy of the matching `#init_states` entry or `default_state`.
Duplicated indicies in a connection can be used to access a specific ID or ID range.

### Time-independent Components

```@docs
TimeIndependentComponent
```

`TimeIndependentComponent(name, func, initial_state)` computes output from an input state rather than advancing with time.
It is useful for algebraic transformations in connector graphs, but does not itself provide a time-evolving step.

It has a single special variable, `#state`.
It is also not possible to define additional variables.

`getstate` returns the `#state`.

## Abstract Types

These abstract types are extension points and are useful for method signatures:

```@docs
AbstractComponent
AbstractTimeDependentComponent
AbstractTimeIndependentComponent
AbstractComponentIntegrator
AbstractMermaidSolver
AbstractMermaidIntegrator
AbstractMermaidProblem
AbstractMermaidSolution
AbstractConnectedVariable
AbstractConnector
```
