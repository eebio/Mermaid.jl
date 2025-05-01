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

    function f1(x, p, t)
        return x - x * p[1]
    end
    function f2(y, p, t)
        return -y + p[1] * y
    end
    prob1 = ODEProblem(f1, u0[1], tspan, [2.0]) # TODO Initial value for params is intentionally wrong
    prob2 = ODEProblem(f2, u0[2], tspan, [4.0])
    c1 = ODEComponent(
        model=prob1,
        name="Prey",
        input_names=["Predator.pop"],
        output_indices=Dict("pop" => 1),
        time_step=0.002,
        alg=Euler(),
        intkwargs=(:adaptive => false,),
    )

    c2 = ODEComponent(
        model=prob2,
        name="Predator",
        time_step=0.002,
        input_names=["Prey.pop"],
        output_indices=Dict("pop" => 1),
        alg=Euler(),
        intkwargs=(:adaptive => false,),
    )

    mp = MermaidProblem(components=[c1, c2], max_t=10.0)

    using CommonSolve # TODO remove the need for this
    alg = MinimumTimeStepper()
    solMer = solve(mp, alg)

    @test all(solMer.t .≈ solODE.t)
    preyODE = [solODE(t)[1] for t in solMer.t]
    predatorODE = [solODE(t)[2] for t in solMer.t]
    @test all(solMer.u["Prey.pop"] .≈ preyODE)
    @test all(solMer.u["Predator.pop"] .≈ predatorODE)
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

    @parameters y
    @variables x(t)
    eqs = [D(x) ~ x - x * y]
    @mtkbuild lv1 = ODESystem(eqs, t)
    prob = ODEProblem(lv1, [x => 4.0], (0.0, 10.0), [y => 2.0])

    c1 = ODEComponent(
        model=prob,
        name="Prey",
        input_names=["Predator.pop"],
        output_indices=Dict("pop" => 1),
        time_step=0.002,
        alg=Euler(),
        intkwargs=(:adaptive => false,),
        mtk_input_symbols=Dict("Predator.pop" => y),
    )

    @parameters x
    @variables y(t)
    eqs = [D(y) ~ -y + x * y]
    @mtkbuild lv2 = ODESystem(eqs, t)
    prob = ODEProblem(lv2, [y => 2.0], (0.0, 10.0), [x => 4.0])

    c2 = ODEComponent(
        model=prob,
        name="Predator",
        time_step=0.002,
        input_names=["Prey.pop"],
        output_indices=Dict("pop" => 1),
        alg=Euler(),
        intkwargs=(:adaptive => false,),
        mtk_input_symbols=Dict("Prey.pop" => x),
    )

    mp = MermaidProblem(components=[c1, c2], max_t=10.0)

    using CommonSolve
    alg = MinimumTimeStepper()
    # Ensure the code is compiled
    solMer = solve(mp, alg)

    @test all(solMer.t .≈ solODE.t)
    preyODE = [solODE(t)[1] for t in solMer.t]
    predatorODE = [solODE(t)[2] for t in solMer.t]
    @test all(solMer.u["Prey.pop"] .≈ preyODE)
    @test all(solMer.u["Predator.pop"] .≈ predatorODE)
end
