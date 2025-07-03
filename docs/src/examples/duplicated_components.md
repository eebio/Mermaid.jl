# Advanced Duplicated Components

You may have seen us use duplicated components in the [Tutorial](@ref).
This is a very powerful tool that lets you efficiently create many instances of a component integrator, each with their own independent state that can be stepped independently.
In the [Tutorial](@ref), we used duplicated components to create lots of instances of the tree ODE model, so every tree could be tracked independently.
Rather than creating 640 components, each with its own integrator - we create 640 state vectors, reducing the memory requirements for the duplicated component.

However, while this functionality is useful, it is not always possible to specify the number of instances *a priori*.
For this reason, it is also possible to create duplicated components with a variable number of instances.

## Setup

For this example, we are going to create an agent-based model of a cell population with each cell goverened by a simple metabolic growth model.

```@example dupcomp
using DifferentialEquations, Agents, Random, Mermaid
function cell!(du, u, p, t)
    nutrients, mass = u
    du[1] = 0
    yield = 1 # TODO
    vmax = 1
    K_nutrients = 1
    decay = 1
    uptake = vmax * nutrients/(K_nutrients + nutrients)
    du[2] = yield*uptake - decay*u[2]
end

u0 = [4.0, 2.0] #TODO
tspan = (0.0, 150.0) #TODO
prob = ODEProblem(cell!, u0, tspan)
using Mermaid
comp1 = ODEComponent(
    model=prob,
    name="cell",
    state_names=Dict("nutrients" => 1, "mass" => 2),
)

@agent struct Cell(ContinuousAgent{2, Float64})
    mass::Float64 # Cell mass is informed by ODE model
    nutrients::Float64 # Local nutrient availability
end

function colony(; n_cells=20, griddims=(40, 40), seed=2)
    space = ContinuousSpace(griddims; periodic=false)
    rng = Random.MersenneTwister(seed)
    colony = StandardABM(Cell, space; rng, agent_step! = cell_step!)
    for _ in 1:n_cells
        vel = rand(2).-0.5
        mass = rand()
        nutrients = rand()
        add_agent!(colony, vel, mass, nutrients)
    end
    return colony
end

function cell_step!(cell, colony)
    # If large mass, split into two
    move_agent!(cell...)
    cell.nutrients = f(cell.pos)
end

pop = colony()

comp2 = AgentsComponent(
    model=pop,
    name="population",
    state_names=Dict("mass" => :mass, "nutrients" => :nutrients)
)

conn1 = Connector(inputs=["colony.nutrients"], outputs=["cell.nutrients"])
conn2 = Connector(inputs=["cell.mass"], outputs=["colony.mass"])
nothing # hide
```

## Flexible Duplicated Components

If we don't provide any value for the `instances` field when creating the duplicated component, it will be created as a flexible duplicated component, letting a special input called `#ids` give the indexes of the current states (any indexes not specified in `#ids` are removed, and any state indexes specified in `#ids` that are not current states are created).
Some components, like AgentComponents, already have pre-made special outputs for the `#ids` which we can use to couple duplicated components to an Agent-based model.

```@example dupcomp
dup_comp = DuplicatedComponent(
    component=comp1,
    initial_states=[copy(u0) for _ in 1:640] # TODO
)
conn3 = Connector(inputs=["colony.#ids"], outputs=["cell.#ids"])

nothing # hide
```

## Solving

```julia
mp = MermaidProblem(components=[dup_comp, comp2], connectors=[conn1, conn2], max_t=tspan[2])
alg = MinimumTimeStepper()
sol = solve(mp, alg)
```
