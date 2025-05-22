@testitem "agent" begin
    using Agents, DifferentialEquations

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
        return 2 * 1 / 8 * cos(t / 8)
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

    @test solMer["ode.happy"][1:end-1] == solMer["Schelling.min_to_be_happy"][2:end]

    intMer = init(mp, alg)
    for _ in 1:10
        step!(intMer, 1.0)
    end
    # Early on, min_to_be_happy is high, so lots of agents moving
    pos = [intMer.integrators[1].integrator[i].pos for i in 1:300]
    step!(intMer, 1.0)
    pos2 = [intMer.integrators[1].integrator[i].pos for i in 1:300]
    @test any(pos != pos2)
    for _ in 1:35
        step!(intMer, 1.0)
    end
    # Later, min_to_be_happy is low, so agents aren't moving
    pos3 = [intMer.integrators[1].integrator[i].pos for i in 1:300]
    step!(intMer, 1.0)
    pos4 = [intMer.integrators[1].integrator[i].pos for i in 1:300]
    @test all(pos3 == pos4)
    # And all agents are happy
    @test all([intMer.integrators[1].integrator[i].mood for i in 1:300])
end
