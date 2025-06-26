using DifferentialEquations
using Mermaid
using CommonSolve
using Statistics

function tree!(du, u, p, t)
    x, y = u
    du[1] = 0
    du[2] = (y * (1 - y / 10.0) - x * y) / 10
end
u0 = [4.0, 2.0]
tspan = (0.0, 150.0)
prob = ODEProblem(tree!, u0, tspan)

comp1 = ODEComponent(
    model=prob,
    name="tree",
    state_names=Dict("heat" => 1, "life" => 2),
)
dup_comp = DuplicatedComponent(
    component=comp1,
    instances=640,
    initial_states=[copy(u0) for _ in 1:640]
)
using Agents, Random
@agent struct Tree(GridAgent{2})
    heat::Float64 # Heat is averaged across neighbors, passed to ODE model
    life::Float64 # Life is informed by ODE model
end

function forest_fire(; density=0.4, griddims=(40, 40), seed=2)
    space = GridSpaceSingle(griddims; periodic=false, metric=:chebyshev)
    rng = Random.MersenneTwister(seed)
    forest = StandardABM(Tree, space; rng, (agent_step!)=tree_step!)
    for _ in 1:floor(density * prod(griddims))
        # Randomly place trees in the grid
        add_agent_single!(forest; heat=rand(), life=rand())
    end
    return forest
end

function tree_step!(tree, forest)
    nearbyheat = mean([getproperty(neighbor, :heat) for neighbor in nearby_agents(tree, forest, 1)])
    if isnan(nearbyheat)
        nearbyheat = 0.0
    end
    tree.heat = tree.heat * 0.9 + nearbyheat * 0.1
    if rand(abmrng(forest)) < 1e-4 # Random chance of fire
        tree.heat = 10.0
    end
    # Simulate tree life cycle
    if tree.heat > 1.0 && tree.life > 1.0
        # Tree on fire
        tree.heat += 1.0
    else
        # Tree not on fire so heat disappates
        tree.heat -= 0.05
    end
    if tree.heat < 0.0
        tree.heat = 0.0
    end
end

forest = forest_fire()

comp2 = AgentsComponent(
    model=forest,
    name="forest",
    state_names=Dict("heat" => :heat, "life" => :life)
)

conn1 = Connector(inputs=["forest.heat[1:640]"], outputs=["tree[1:640].heat"])
conn2 = Connector(inputs=["tree[1:640].life"], outputs=["forest.life[1:640]"])
using Makie
using CairoMakie

groupcolor(tree) = tree.heat > 1 ? :red : :green
groupmarker(a) = a.life > 1 ? :utriangle : :circle
fig, ax = abmplot(forest; agent_color=groupcolor, agent_marker=groupmarker, agent_size=10)
io = VideoStream(fig)
function plot_input(model)
    empty!(ax)
    abmplot!(ax, model; agent_color=groupcolor, agent_marker=groupmarker, agent_size=10)
    recordframe!(io)
end

conn3 = Connector(
    inputs=["forest.#model"],
    outputs=Vector{String}(),
    func=(model) -> plot_input(model)
)

mp = MermaidProblem(components=[dup_comp, comp2], connectors=[conn1, conn2, conn3], max_t=tspan[2])
alg = MinimumTimeStepper()
sol = solve(mp, alg)

save("examples/outputs/forest_fire.mp4", io)
