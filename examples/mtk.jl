using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
using Mermaid

maxt = [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]
odeTimes = []
mermaidTimes = []

for tend in maxt
    # Using DifferentialEquations.jl for an ODE system
    # Lotka-Volterra equations
    @variables x(t) y(t)
    eqs = [D(x) ~ x - x * y
        D(y) ~ -y + x * y]
    @mtkbuild lv = ODESystem(eqs, t)
    prob = ODEProblem(lv, [x => 4.0, y => 2.0], (0.0, tend))

    global solODE = solve(prob, Euler(); adaptive=false, dt=0.002)
    plot(solODE, title="Lotka-Volterra ODE", xlabel="Time", ylabel="Population", label=["Prey ODE" "Predator ODE"])
    push!(odeTimes, @elapsed solve(prob, Euler(); adaptive=false, dt=0.002))

    eqs = [D(x) ~ x - x * y
        D(y) ~ 0]
    @mtkbuild lv1 = ODESystem(eqs, t)
    prob = ODEProblem(lv1, [x => 4.0, y => 2.0], (0.0, tend))

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
    prob = ODEProblem(lv2, [x => 4.0, y => 2.0], (0.0, tend))

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

    mp = MermaidProblem(components=[c1, c2], connectors=[conn1, conn2], max_t=tend)

    using CommonSolve
    alg = MinimumTimeStepper()
    # Ensure the code is compiled
    global solMer = solve(mp, alg)
    plot!(solMer.t, solMer.u["Prey.prey"], label="Prey Mermaid")
    display(plot!(solMer.t, solMer.u["Predator.predator"], label="Predator Mermaid"))
    push!(mermaidTimes, @elapsed solve(mp, alg))
end

@show mermaidTimes ./ odeTimes
