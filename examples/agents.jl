using Agents, DifferentialEquations, Mermaid

space = GridSpace((20, 20))

@agent struct Schelling(GridAgent{2})
    mood::Bool = false
    group::Int
end

function schelling_step!(agent, model)
    minhappy = model.min_to_be_happy
    count_neighbors_same_group = 0
    for neighbor in nearby_agents(agent, model)
        if agent.group == neighbor.group
            count_neighbors_same_group += 1
        end
    end
    if count_neighbors_same_group ≥ minhappy
        agent.mood = true
    else
        agent.mood = false
        move_agent_single!(agent, model)
    end
    return
end

properties = Dict(:min_to_be_happy => 3.0)

model = StandardABM(
    Schelling,
    space;
    (agent_step!)=schelling_step!, properties
)

for n in 1:300
    add_agent_single!(model; group=n < 300 / 2 ? 1 : 2)
end

using Mermaid

c1 = AgentsComponent(
    model=model,
    name="Schelling",
    state_names=Dict("min_to_be_happy" => :min_to_be_happy),
    time_step=1.0,
)

function f2(u, p, t)
    return 2*1/8*cos(t/8)
end
u0 = 3.0
tspan = (0.0, 100.0)
prob = ODEProblem(f2, u0, tspan)
c2 = ODEComponent(
    model=prob,
    name="ode",
    time_step=1.0,
    state_names=Dict("happy" => 1),
    alg=Tsit5(),
)

conn = Connector(
    inputs=["ode.happy"],
    outputs=["Schelling.min_to_be_happy"],
)

mp = MermaidProblem(components=[c1, c2], connectors=[conn], max_t=100.0)

using CommonSolve
alg = MinimumTimeStepper()
intMer = init(mp, alg)
solMer = solve!(intMer)

plot(solMer.t, solMer.u["ode.happy"])

using CairoMakie

groupcolor(a) = a.mood ? :blue : :orange
groupmarker(a) = a.group == 1 ? :circle : :rect

intMer = init(mp, alg)
for i in 1:100
    figure, ax = abmplot(intMer.integrators[1].integrator; agent_color=groupcolor, agent_marker=groupmarker, agent_size=10)
    ax.title = "Min to be happy: $(intMer.integrators[1].integrator.min_to_be_happy)"
    display(figure) # returning the figure displays it
    step!(intMer, 1.0)
end
