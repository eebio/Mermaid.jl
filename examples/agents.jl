using Agents

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

properties = Dict(:min_to_be_happy => 3)

model = StandardABM(
    Schelling,
    space;
    (agent_step!)=schelling_step!, properties
)

for n in 1:300
    add_agent_single!(model; group=n < 300 / 2 ? 1 : 2)
end

using Mermaid

comp = AgentsComponent(
    model=model,
    name="Schelling",
    state_names=Dict("min_to_be_happy" => :min_to_be_happy),
    time_step=1.0,
)

mp = MermaidProblem(components=[comp], connectors=[], max_t=100.0)

using CommonSolve
alg = MinimumTimeStepper()
# Ensure the code is compiled
global solMer = solve(mp, alg)
