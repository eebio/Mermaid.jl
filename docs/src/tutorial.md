# Tutorial

In this tutorial, we will:

* Create a hybrid simulation between an Agent-based model, defined in Agents.jl, and an ODE system defined through DifferentialEquations.jl.
* Introduce Mermaid Components for Agents.jl and DifferentialEquations.jl.
* Demonstrate how these Components can be connected together through Connections.
* Solve the hybrid model.
* Visualise the results of the simulation.

The example system we will use for this will be a model of a forest fire (governed by an Agent-based model), with the growth of each tree informed by an ODE model.

## Components

To begin, we need to define our components. These will be an ODE model component for each tree, and an Agent-based model component for handling the forest level properties (in this case, the spread of heat/fire).

### ODE Components

To define the ODE model, let's have a look at how to define an ODE Component.

```@docs; canonical=false
ODEComponent
```

We see that we need to define an `ODEProblem` to use in the component, so let's create one.

```@example tutorial
using DifferentialEquations
function tree!(du, u, p, t)
    x, y = u
    du[1] = x - x * y
    du[2] = 0
end
u0 = [4.0,2.0]
tspan = (0.0, 10.0)
prob = ODEProblem(tree!, u0, tspan)
```

Next, we want to wrap this ODEProblem inside an ODEComponent.
For this, we will need to define the `state_names` field, and should generally provide a value for the `name` field (since component names in a hybrid simulation should be unique).

```@example tutorial
using Mermaid
comp1 = ODEComponent(
    model=prob,
    name="tree",
    state_names=Dict("heat" => 1, "life" => 2),
)
```

### Duplicated Components

Since the `ODEProblem` is defined for only a single tree, we can efficiently simulate the ODE by generating a duplicated component.
This component stores a single `ODEProblem` will solve across many different states.
In this case, we can have a state for each tree in the Agent-based model.
Let's have a look at how to define a [DuplicatedComponent](@ref).

```@docs; canonical=false
DuplicatedComponent
```

```julia
dup_comp = DuplicatedComponent(
    component=comp1,
    instances=640,
    states=repeat([u0], 640)
)
```

### Agents.jl Components

Now that we have created our [ODEComponent](@ref), we can move on to the [AgentsComponent](@ref), so let's have a look at its documentation.

```@docs; canonical=false
AgentsComponent
```

We can see that, again, we need to define the model (this time a `StandardABM` from `Agents.jl`), a `name` and a `state_names`.

```@example tutorial
using Agents, Random
@agent struct Tree(GridAgent{2})
    heat::Float64 # Heat is averaged across neighbors, passed to ODE model
    life::Float64 # Life is informed by ODE model
end

function forest_fire(; density=0.4, griddims=(40, 40), seed=2)
    space = GridSpaceSingle(griddims; periodic=false, metric=:chebyshev)
    rng = Random.MersenneTwister(seed)
    forest = StandardABM(Tree, space; rng, agent_step! = tree_step!)
    for _ in 1:floor(density * prod(griddims))
        # Randomly place trees in the grid
        add_agent_single!(forest; heat=rand(), life=rand())
    end
    return forest
end

function tree_step!(tree, forest)
    tree.heat = mean([getproperty(neighbor, :heat) for neighbor in nearby_agents(tree, forest, 1)])
    if rand(abmrng(forest)) < 1e-4 # Random chance of fire
        tree.heat = 10.0
    end
    # Simulate tree life cycle
    if tree.heat > 1.0 && tree.life > 0
        # Tree on fire
        tree.heat += 1.0
    else
        # Tree not on fire so heat disappates
        tree.heat -= 0.1
    end
end

forest = forest_fire()

comp2 = AgentsComponent(
    model=forest,
    name="forest",
    state_names=Dict("heat" => :heat, "life" => :life)
)
```

## Connections

We can now set up the connections between the variables in the two components.

```@docs; canonical=false
Connector
```

```julia
conn1 = Connection(inputs=["forest.heat[1:640]"], outputs=["tree[1:640].heat"])
conn2 = Connection(inputs=["tree[i]."], outputs=[""])
```

## Solving the hybrid model

To create the hybrid model, we need to create a [MermaidProblem](@ref).
We can then solve this using the `CommonSolve` interface.

```@docs; canonical=false
MermaidProblem
```

```julia
mp = MermaidProblem(components = [comp1, dup_comp2], connectors = [conn1, conn2])
alg = MinimumTimeStepper()
sol = solve(mp, alg)
```

## Plotting the solution

After running `solve`, we get `sol`, a [MermaidSolution](@ref) instance.
This stores all variables given in `state_names` at each timepoint.

```@docs; canonical=false
MermaidSolution
```

```julia
using Plots

plot(sol.t, sol["forest.life[1]"])
```
