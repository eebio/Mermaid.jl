@testitem "simple ODE" begin
    using DifferentialEquations

    function f!(du, u, p, t)
        x, y = u
        du[1] = x - x * y
        du[2] = -y + x * y
    end

    u0 = [4.0, 2.0]
    tspan = (0.0, 10.0)
    prob = ODEProblem(f!, u0, tspan)
    solODE = solve(prob, Euler(); adaptive=false, dt=0.002)

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
    prob1 = ODEProblem(f1!, [u0[1],2.0], tspan) # TODO Initial value for params is intentionally wrong
    prob2 = ODEProblem(f2!, [u0[2],4.0], tspan)
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

    using CommonSolve # TODO remove the need for this
    alg = MinimumTimeStepper()
    solMer = solve(mp, alg)

    @test all(solMer.t .≈ solODE.t)
    preyODE = [solODE(t)[1] for t in solMer.t]
    predatorODE = [solODE(t)[2] for t in solMer.t]
    @test all(solMer.u["Prey.prey"] .≈ preyODE)
    @test all(solMer.u["Predator.predator"] .≈ predatorODE)
end

@testitem "mtk" begin
    using ModelingToolkit, DifferentialEquations
    using ModelingToolkit: t_nounits as t, D_nounits as D

    @variables x(t) y(t)
    eqs = [D(x) ~ x - x * y
        D(y) ~ -y + x * y]
    @mtkbuild lv = ODESystem(eqs, t)
    prob = ODEProblem(lv, [x => 4.0, y => 2.0], (0.0, 10.0))

    solODE = solve(prob, Euler(); adaptive=false, dt=0.002)

    eqs = [D(x) ~ x - x * y
        D(y) ~ 0]
    @mtkbuild lv1 = ODESystem(eqs, t)
    prob = ODEProblem(lv1, [x => 4.0, y => 2.0], (0.0, 10.0))

    c1 = ODEComponent(
        model=prob,
        name="Prey",
        time_step=0.002,
        state_names=Dict("prey" => x, "predator" => y),
        alg=Euler(),
        intkwargs=(:adaptive => false,),
    )

    eqs = [D(x) ~ 0
        D(y) ~ -y + x * y]
    @mtkbuild lv2 = ODESystem(eqs, t)
    prob = ODEProblem(lv2, [x => 4.0, y => 2.0], (0.0, 10.0))

    c2 = ODEComponent(
        model=prob,
        name="Predator",
        time_step=0.002,
        state_names=Dict("prey" => x, "predator" => y),
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

    using CommonSolve
    alg = MinimumTimeStepper()
    # Ensure the code is compiled
    solMer = solve(mp, alg)

    @test all(solMer.t .≈ solODE.t)
    preyODE = [solODE(t)[1] for t in solMer.t]
    predatorODE = [solODE(t)[2] for t in solMer.t]
    @test all(solMer.u["Prey.prey"] .≈ preyODE)
    @test all(solMer.u["Predator.predator"] .≈ predatorODE)
end
