using DifferentialEquations, Agents, Random, Mermaid, CairoMakie
using LinearAlgebra: norm

Random.seed!(1)

function cell!(du, u, p, t)
    nutrients = u[1]
    du[1] = 0
    uptake = 2 * nutrients / (1 + nutrients)
    decay = 0.1 * (1+u[2]/10)
    du[2] = uptake - decay
end

u0 = [1.0, 1.0]
tspan = (0.0, 250.0)
prob = ODEProblem(cell!, u0, tspan)
using Mermaid
comp1 = ODEComponent(
    model=prob,
    name="cell",
    state_names=Dict("nutrients" => 1, "mass" => 2),
)

@agent struct Cell(ContinuousAgent{2,Float64})
    mass::Float64 # Cell mass is informed by ODE model
    nutrients::Float64 # Local nutrient availability
end

function colony(; n_cells=3, n_nodes=5, griddims=(40, 40), seed=2)
    space = ContinuousSpace(griddims; periodic=true)
    rng = Random.MersenneTwister(seed)
    nodes = [rand(2) .* griddims for _ in 1:n_nodes]
    colony = StandardABM(Cell, space; rng, (agent_step!)=cell_step!, properties=Dict(:nodes => nodes))
    for _ in 1:n_cells
        vel = rand(2) .- 0.5
        mass = rand()
        nutrients = rand()
        add_agent!(colony, vel, 1, 1)
    end
    return colony
end

function wrap_periodic(pos, dims)
    return SVector{length(pos)}(mod.(pos, dims))
end

function periodic_distance(a, b, dims)
    # Computes minimum image distance between points a and b in periodic box of size dims
    d = abs.(a .- b)
    return norm(min.(d, dims .- d))
end

function nutrients(pos, colony)
    nutrients = 0.0
    spread = 7.0
    dims = abmspace(colony).extent
    for node in colony.nodes
        # Calculate periodic distance to node
        dist = periodic_distance(node, pos, dims)
        nutrients += exp(-dist^2 / spread^2)/2 # Gaussian decay
    end
    return nutrients
end

function cell_step!(cell, colony)
    # Move away from nearby agents
    speed = 0.5
    for cell2 in nearby_agents(cell, colony, 0.5)
        cell.vel -= speed * (cell2.pos - cell.pos) / norm(cell2.pos - cell.pos)^2
    end
    # Max speed
    if norm(cell.vel) > speed
        cell.vel = speed * cell.vel / norm(cell.vel)
    end
    # Walk and apply chemotaxis
    oldnutrients = nutrients(cell.pos, colony)
    walk!(cell, cell.vel, colony)
    newnutrients = nutrients(cell.pos, colony)
    if newnutrients < oldnutrients
        # If nutrients decrease, random direction on next iteration
        cell.vel += rand(2) .- 0.5
    end
    # Update nutrients and mass
    cell.nutrients = nutrients(cell.pos, colony)/(length(collect(nearby_ids(cell, colony, 0.5)))+1)
    # If large mass, split into two
    splitmass = 15
    if cell.mass > splitmass
        dims = abmspace(colony).extent
        for m in (splitmass/2, cell.mass - splitmass/2)
            new_pos = wrap_periodic(cell.pos + rand(2) .- 0.5, dims)
            add_agent!(new_pos, colony; vel=cell.vel, mass=m, nutrients=0)
        end
        remove_agent!(cell, colony)
    end
    if cell.mass < 0
        remove_agent!(cell, colony)
    end
end

function make_nutrient_heatarray(colony; gridsize=(40, 40))
    arr = zeros(Float64, gridsize...)
    xs = range(0, stop=abmspace(colony).extent[1], length=gridsize[1])
    ys = range(0, stop=abmspace(colony).extent[2], length=gridsize[2])
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
        arr[i, j] = nutrients(SVector(x, y), colony)
    end
    return arr
end

pop = colony()

comp2 = AgentsComponent(
    model=pop,
    name="colony",
    state_names=Dict("mass" => :mass, "nutrients" => :nutrients)
)

conn1 = Connector(inputs=["colony.nutrients"], outputs=["cell.nutrients"])
conn2 = Connector(inputs=["cell.mass"], outputs=["colony.mass"])

dup_comp = DuplicatedComponent(
    component=comp1,
    initial_states=[],
    default_state=u0,
)
conn3 = Connector(inputs=["colony.#ids"], outputs=["cell.#ids"])

fig, ax = abmplot(pop; agent_color=:black, agent_marker=:circle, agent_size=x->x.mass+3,
    heatarray=make_nutrient_heatarray, heatkwargs=(colormap=:Greens_5,))
io = VideoStream(fig; framerate=10)
function plot_input(model)
    empty!(ax)
    abmplot!(ax, model; agent_color=:black, agent_marker=:circle, agent_size=x -> x.mass + 3,
        heatarray=make_nutrient_heatarray, heatkwargs=(colormap=:Greens_5,))
    ax.title = "Time: $(round(abmtime(model))), Population: $(nagents(model))"
    recordframe!(io)
end

conn4 = Connector(
    inputs=["colony.#model"],
    outputs=Vector{String}(),
    func=(model) -> plot_input(model)
)

mp = MermaidProblem(components=[dup_comp, comp2], connectors=[conn1, conn2, conn3, conn4], max_t=tspan[2])
alg = MinimumTimeStepper()
sol = solve(mp, alg)

save("examples/outputs/cell_colony.mp4", io)
