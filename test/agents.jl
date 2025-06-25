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

@testitem "state control" begin
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

    properties = Dict(:min_to_be_happy => 3.0, :list_property => [1, 2, 3])

    model = StandardABM(
        Schelling,
        space;
        (agent_step!)=schelling_step!, properties
    )

    for n in 1:300
        add_agent_single!(model; group=n < 300 / 2 ? 1 : 2)
    end

    c1 = AgentsComponent(
        model=model,
        name="Schelling",
        state_names=Dict("min_to_be_happy" => :min_to_be_happy, "list_property" => :list_property,"mood" => :mood, "group" => :group),
        time_step=0.2,
    )

    conn1 = Connector(
        inputs=["Schelling.min_to_be_happy"],
        outputs=["other.min_to_be_happy"],
    )
    conn2 = Connector(
        inputs=["other.min_to_be_happy"],
        outputs=["Schelling.group[1:300]"],
    )

    integrator = init(c1, [conn1, conn2])

    @test issetequal(Mermaid.variables(integrator), ["min_to_be_happy", "group", "mood", "list_property"])

    # Check initial state
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.min_to_be_happy")) == 3.0
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.group")) == [n < 300 / 2 ? 1 : 2 for n in allids(integrator.integrator)]
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.mood")) == [false for _ in allids(integrator.integrator)]
    # Check setting state
    Mermaid.setstate!(integrator, ConnectedVariable("Schelling.min_to_be_happy"), 5.0)
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.min_to_be_happy")) == 5.0
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.group")) == [n < 300 / 2 ? 1 : 2 for n in allids(integrator.integrator)]
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.mood")) == [false for _ in allids(integrator.integrator)]

    Mermaid.setstate!(integrator, ConnectedVariable("Schelling.group[1:3]"), [3, 4, 5])
    Mermaid.setstate!(integrator, ConnectedVariable("Schelling.mood[2:3]"), [true, true])
    Mermaid.setstate!(integrator, ConnectedVariable("Schelling.mood[300]"), true)
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.min_to_be_happy")) == 5.0
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.group")) == [n ∈ [1,2,3] ? n+2 : (n < 300 / 2 ? 1 : 2) for n in allids(integrator.integrator)]
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.mood")) == [n ∈ [2,3, 300] ? true : false for n in allids(integrator.integrator)]
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.mood[300]")) == true

    # Check time control (can't set time in Agents.jl)
    @test Mermaid.gettime(integrator) == 0.0
    step!(integrator)
    @test Mermaid.gettime(integrator) == 0.2

    # Step means the state has changed
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.group")) ≠ [3, 4, 5, [n < 300 / 2 ? 1 : 2 for n in 4:300]...]
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.mood")) ≠ [false, true, true, [false for _ in 4:300]..., true]

    Mermaid.setstate!(integrator, ConnectedVariable("Schelling.mood"), [true for _ in allids(integrator.integrator)])
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.mood")) == [true for _ in allids(integrator.integrator)]
    Mermaid.setstate!(integrator, ConnectedVariable("Schelling.list_property"), [4, 5, 6])
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.list_property")) == [4, 5, 6]
    Mermaid.setstate!(integrator, ConnectedVariable("Schelling.list_property[1:2]"), [7, 8])
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.list_property[1]")) == 7
    @test Mermaid.getstate(integrator, ConnectedVariable("Schelling.list_property[1:3]")) == [7, 8, 6]

    # getstate and setstate! for duplicated AgentsComponent
    @test Mermaid.getstate(integrator) isa StandardABM
    @test Mermaid.gettime(integrator) == 0.2
    state = Mermaid.getstate(integrator, true) # Get a copy of the state
    state2 = Mermaid.getstate(integrator) # Default is don't copy, just return reference
    step!(integrator)
    @test Mermaid.gettime(integrator) == 0.4
    Mermaid.setstate!(integrator, state)
    @test Mermaid.gettime(integrator) == 0.2
    Mermaid.setstate!(integrator, state2)
    @test Mermaid.gettime(integrator) == 0.4

    # Settime does nothing since Agents.jl does not support setting time directly but rather stores it within the state
    Mermaid.settime!(integrator, 1.0)
    @test Mermaid.gettime(integrator) == 0.4
end
