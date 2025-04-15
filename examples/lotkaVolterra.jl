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
label=["Prey" "Predator"], legend=:topright)

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
    label=["Prey" "Predator"], legend=:topright, linestyle=:dot)
