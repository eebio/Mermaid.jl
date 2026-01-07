using OrdinaryDiffEq, Random, Mermaid, BenchmarkTools, Plots, LinearAlgebra, Statistics
using Surrogates, Flux
# Random.seed!(1)

using Profile

function lorenz!(du, u, p, t)
    freq = 1.0 #1000.0 makes surrogates 15x faster
    du[1] = sin(freq * t*2pi) + cos(freq * t*2pi)
    du[2] = cos(freq * t*2pi)
    du[3] = sin(freq * t*2pi)
    nothing
end

u0 = [1.0, 1.0, 0.5]
tspan = (0.0, 250.0)
prob = ODEProblem(lorenz!, u0, tspan)

comp1 = DEComponent(prob, Tsit5();
    name = "cell", time_step = 1.0, state_names = Dict("nutrients" => 1, "mass" => 2, "other" => 3),
    intkwargs = (abstol = 1e-10, reltol = 1e-10, save_everystep = false)
)

# Create a surrogate component
surrogate_comp = SurrogateComponent(comp1, [0.4, 0.7, 0.1], [1.3, 1.5, 1.1];
    n_samples = 10000, n_epochs=10000)

println("=== Timing Comparison ===")

# Set up initial conditions for timing test
alg = MinimumTimeStepper()

# Benchmark the original component
println("Benchmarking original component...")
mp1 = MermaidProblem(components = [comp1], connectors = Connector[], max_t = 100.0)
sol1 = solve(mp1, alg) # Compile first

original_time = @benchmark solve!(int) setup=(int = init(
    $mp1, $alg; save_vars = ["cell.nutrients"]))

# Benchmark the surrogate component
println("Benchmarking surrogate component...")
mp2 = MermaidProblem(components = [surrogate_comp], connectors = Connector[], max_t = 100.0)
sol2 = solve(mp2, alg) # Compile first

surrogate_time = @benchmark solve!(int) setup=(int = init(
    $mp2, $alg; save_vars = ["cell.nutrients"]))

println("Original component minimum time: ", minimum(original_time.times) / 1e6, " ms")
println("Surrogate component minimum time: ", minimum(surrogate_time.times) / 1e6, " ms")
println("Speedup: ", minimum(original_time.times) / minimum(surrogate_time.times), "x")

println("\n=== Visualization ===")
# Plot comparison
p1 = plot(sol1.t, sol1["cell.nutrients"], label = "Original - Nutrients", linewidth = 2)
plot!(p1, sol2.t, sol2["cell.nutrients"], label = "Surrogate - Nutrients", linestyle = :dash, linewidth = 2)
plot!(p1, sol1.t, sol1["cell.mass"], label = "Original - Mass", linewidth = 2)
plot!(p1, sol2.t, sol2["cell.mass"], label = "Surrogate - Mass", linestyle = :dash, linewidth = 2)
plot!(p1, sol1.t, sol1["cell.other"], label = "Original - Other", linewidth = 2)
plot!(p1, sol2.t, sol2["cell.other"], label = "Surrogate - Other", linestyle = :dash, linewidth = 2)
xlabel!(p1, "Time")
