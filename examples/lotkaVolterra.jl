using DifferentialEquations
using Plots

# Using DifferentialEquations.jl for an ODE system
# Lotka-Volterra equations
function f!(du, u, p, t)
    x, y = u
    du[1] = x - x * y
    du[2] = -y + x * y
end

# Set up ODE Problem
u0 = [4.0, 2.0]
tspan = (0.0, 10.0)
prob = ODEProblem(f!, u0, tspan)
@time sol = solve(prob, Euler(); adaptive=false, dt=0.01)

# Plot the solution
plot(sol, title="Lotka-Volterra Model", xlabel="Time", ylabel="Population",
label=["Prey 1" "Predator 1"], legend=:topright)

# Create the system as components
using Mermaid

# Define the ODE component
# Variables that are held constant should be defined as parameters
# Prey
function f1(x, p, t)
    return x - x * p[1]
end
x0 = 4.0
tspan = (0.0, 10.0)
prob = ODEProblem(f1, x0, tspan)
c1 = ODEComponent(
    model=prob,
    name="Prey",
    outputs=Dict("pop" => x0),
    inputs=Dict("Predator.pop" => 0.0),
    output_indices=Dict("pop" => 1),
    state=x0,
    time_step=0.01,
    )

# Predator
function f2(y, p, t)
    return -y + p[1] * y
end
y0 = 2.0
tspan = (0.0, 10.0)
prob = ODEProblem(f2, y0, tspan)
c2 = ODEComponent(model=prob, name="Predator", time_step=0.01, state=y0, output_indices=Dict("pop" => 1),
    inputs=Dict("Prey.pop" => 0.0), outputs=Dict("pop" => y0))

# Create the system of components
components = Vector{ODEComponent}([c1, c2])

@time components, outputs = Mermaid.solve!(components, 10.0)

plot!(outputs[:,1], outputs[:,2:3], title="Lotka-Volterra Model", xlabel="Time", ylabel="Population",
    label=["Prey 2" "Predator 2"], legend=:topright, linestyle=:dot)



# CommonSolve Mermaid

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
    time_step=0.01,
)

# Predator
function f2(y, p, t)
    return -y + p[1] * y
end
y0 = 2.0
tspan = (0.0, 10.0)
prob = ODEProblem(f2, y0, tspan, [4.0])
c2 = ODEComponent(model=prob, name="Predator", time_step=0.01, state=y0, output_indices=Dict("pop" => 1),
    inputs=Dict("Prey.pop" => 0.0), outputs=Dict("pop" => y0))

mp = MermaidProblem(components=[c1,c2], max_t=10.0)

using CommonSolve
alg = MermaidSolver()
@time sol = solve(mp, alg)
plot!(sol.t, sol.u["Prey.pop"], label="Prey 3", linestyle=:dash)
plot!(sol.t, sol.u["Predator.pop"], label="Predator 3", linestyle=:dash)
# TODO we should make sure we stop at the maxt too
