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
sol = solve(prob, Tsit5(); adaptive=false, dt=0.1)

# Plot the solution
plot(sol, title="Lotka-Volterra Model", xlabel="Time", ylabel="Population",
label=["Prey" "Predator"], legend=:topright)

# Create the system as components
