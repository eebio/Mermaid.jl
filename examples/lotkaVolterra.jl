using DifferentialEquations
using Plots
using Mermaid

maxt = [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]
odeTimes = []
mermaidTimes = []

for tend in maxt
    # Using DifferentialEquations.jl for an ODE system
    # Lotka-Volterra equations
    function f!(du, u, p, t)
        x, y = u
        du[1] = x - x * y
        du[2] = -y + x * y
    end

    # Set up ODE Problem
    u0 = [4.0, 2.0]
    tspan = (0.0, tend)
    prob = ODEProblem(f!, u0, tspan)
    global solODE = solve(prob, Euler(); adaptive=false, dt=0.002)
    plot(solODE, title="Lotka-Volterra ODE", xlabel="Time", ylabel="Population", label=["Prey ODE" "Predator ODE"])
    push!(odeTimes, @elapsed solve(prob, Euler(); adaptive=false, dt=0.002))

    function f1!(du, u, p, t)
        x, y = u
        du[1] = x - x * y
        du[2] = 0
    end
    u0 = [4.0,2.0]
    tspan = (0.0, 10.0)
    prob = ODEProblem(f1!, u0, tspan) # Initial value for params is intentionally wrong
    c1 = ODEComponent(
        model=prob,
        name="Prey",
        time_step=0.002,
        state_names=Dict("prey" => 1, "predator" => 2),
        alg=Euler(),
        intkwargs=(:adaptive => false,),
    )

    # Predator
    function f2!(du, u, p, t)
        y, x = u
        du[1] = -y + x * y
        du[2] = 0
    end
    u0 = [2.0,4.0]
    tspan = (0.0, 10.0)
    prob = ODEProblem(f2!, u0, tspan)
    c2 = ODEComponent(
        model=prob,
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

    mp = MermaidProblem(components=[c1, c2], connectors=[conn1, conn2], max_t=tend)

    using CommonSolve
    alg = MinimumTimeStepper()
    # Ensure the code is compiled
    global solMer = solve(mp, alg)
    plot!(solMer.t, solMer.u["Prey.prey"], label = "Prey Mermaid")
    display(plot!(solMer.t, solMer.u["Predator.predator"], label = "Predator Mermaid"))
    push!(mermaidTimes, @elapsed solve(mp, alg))
end

@show mermaidTimes./odeTimes
