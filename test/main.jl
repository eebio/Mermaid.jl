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
