using DifferentialEquations
using Plots

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
    sol = solve(prob, Euler(); adaptive=false, dt=0.002)
    push!(odeTimes, @elapsed sol = solve(prob, Euler(); adaptive=false, dt=0.002))

    # Create the system as components
    using Mermaid

    function f1(x, p, t)
        return x - x * p[1]
    end
    x0 = 4.0
    tspan = (0.0, 10.0)
    prob = ODEProblem(f1, x0, tspan, [2.0]) # Initial value for params is intentionally wrong
    c1 = ODEComponent(
        model=prob,
        name="Prey",
        outputs=Dict("pop" => x0),
        inputs=Dict("Predator.pop" => 0.0),
        output_indices=Dict("pop" => 1),
        state=x0,
        time_step=0.002,
    )

    # Predator
    function f2(y, p, t)
        return -y + p[1] * y
    end
    y0 = 2.0
    tspan = (0.0, 10.0)
    prob = ODEProblem(f2, y0, tspan, [4.0])
    c2 = ODEComponent(model=prob, name="Predator", time_step=0.002, state=y0, output_indices=Dict("pop" => 1),
        inputs=Dict("Prey.pop" => 0.0), outputs=Dict("pop" => y0))

    mp = MermaidProblem(components=[c1,c2], max_t=tend)

    using CommonSolve
    alg = MermaidSolver()
    # Ensure the code is compiled
    sol = solve(mp, alg)
    push!(mermaidTimes, @elapsed sol = solve(mp, alg))
end

@show mermaidTimes./odeTimes
