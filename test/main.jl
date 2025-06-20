@testitem "fullname" begin
    @test Mermaid.fullname(ConnectedVariable("comp", "var", 1:5, [1, 3, 5, 7])) == "comp[[1, 3, 5, 7]].var[1:5]"
    @test Mermaid.fullname(ConnectedVariable("cmp", "var", 1:5, nothing)) == "cmp.var[1:5]"
    @test Mermaid.fullname(ConnectedVariable("cp", "vr", nothing, [1, 3, 5, 7])) == "cp[[1, 3, 5, 7]].vr"
    @test Mermaid.fullname(ConnectedVariable("comp", "var", nothing, nothing)) == "comp.var"
end

@testitem "solution" begin
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
        state_names=Dict("min_to_be_happy" => :min_to_be_happy, "list_property" => :list_property, "mood" => :mood, "group" => :group),
        time_step=1.0,
    )

    mp = MermaidProblem(components=[c1], connectors=[], max_t=10.0)

    using CommonSolve
    alg = MinimumTimeStepper()
    sol = solve(mp, alg)
    @test sol["Schelling.min_to_be_happy"] == [3.0 for _ in sol.t]
    @test sol["Schelling.list_property"] == [[1, 2, 3] for _ in sol.t]
    @test sol["Schelling.list_property[2:3]"] == [[2, 3] for _ in sol.t]
end

@testitem "mermaid integrator" begin
    using DifferentialEquations

    function f1!(du, u, p, t)
        x, y = u
        du[1] = x - x * y
        du[2] = 0
    end
    function f2!(du, u, p, t)
        y, x = u
        du[1] = -y + x * y
        du[2] = 0
    end
    u0 = [4.0, 2.0]
    tspan = (0.0, 10.0)
    prob1 = ODEProblem(f1!, [u0[1], 2.0], tspan) # TODO Initial value for params is intentionally wrong
    prob2 = ODEProblem(f2!, [u0[2], 4.0], tspan)
    c1 = ODEComponent(
        model=prob1,
        name="Prey",
        time_step=0.002,
        state_names=Dict("prey" => 1, "predator" => 2),
        alg=Euler(),
        intkwargs=(:adaptive => false,),
    )

    c2 = ODEComponent(
        model=prob2,
        name="Predator",
        time_step=0.002,
        state_names=Dict("predator" => 1, "prey" => 2),
        alg=Euler(),
        intkwargs=(:adaptive => false,),
    )

    conn1 = Connector(
        inputs=["Predator.predator"],
        outputs=["Prey.predator"],
    )
    conn2 = Connector(
        inputs=["Prey.prey"],
        outputs=["Predator.prey"],
    )

    mp = MermaidProblem(components=[c1, c2], connectors=[conn1, conn2], max_t=10.0)
    integrator = init(mp, MinimumTimeStepper())

    # State control
    @test Mermaid.getstate(integrator, ConnectedVariable("Prey.prey")) == 4.0
    @test Mermaid.getstate(integrator, ConnectedVariable("Predator.predator")) == 2.0
    Mermaid.setstate!(integrator, ConnectedVariable("Prey.prey"), 5.0)
    @test Mermaid.getstate(integrator, ConnectedVariable("Prey.prey")) == 5.0
    @test Mermaid.getstate(integrator, ConnectedVariable("Predator.predator")) == 2.0
    @test Mermaid.getstate(integrator.integrators[1], ConnectedVariable("Prey.prey")) == 5.0
    step!(integrator, 0.01)
    @test Mermaid.getstate(integrator, ConnectedVariable("Prey.prey")) ≠ 5.0
    @test Mermaid.getstate(integrator, ConnectedVariable("Predator.predator")) ≠ 2.0

    # update_inputs!
    integrator = init(mp, MinimumTimeStepper())
    Mermaid.update_inputs!(integrator)
    @test integrator.integrators[1].inputs[ConnectedVariable("Prey.predator")] == 2.0
    @test integrator.integrators[2].inputs[ConnectedVariable("Predator.prey")] == 4.0
    conn1 = Connector(
        inputs=["Predator.predator"],
        outputs=["Prey.predator"],
        func=x->x*4,
    )
    conn2 = Connector(
        inputs=["Prey.prey"],
        outputs=["Predator.prey"],
        func=x->x/1.5,
    )

    mp = MermaidProblem(components=[c1, c2], connectors=[conn1, conn2], max_t=10.0)
    integrator = init(mp, MinimumTimeStepper())
    Mermaid.update_inputs!(integrator)
    @test integrator.integrators[1].inputs[ConnectedVariable("Prey.predator")] == 8.0
    @test integrator.integrators[2].inputs[ConnectedVariable("Predator.prey")] == 4.0/1.5
    conn1 = Connector(
        inputs=["Predator.predator", "Predator.prey"],
        outputs=["Prey.predator", "Prey.prey"],
        func=(x,y) -> x*y,
    )

    mp = MermaidProblem(components=[c1, c2], connectors=[conn1], max_t=10.0)
    integrator = init(mp, MinimumTimeStepper())
    Mermaid.setstate!(integrator, ConnectedVariable("Predator.predator"), 2.0)
    Mermaid.setstate!(integrator, ConnectedVariable("Predator.prey"), 4.0)
    Mermaid.update_inputs!(integrator)
    @test integrator.integrators[1].inputs[ConnectedVariable("Prey.predator")] == 8.0
    @test integrator.integrators[1].inputs[ConnectedVariable("Prey.prey")] == 8.0
end
