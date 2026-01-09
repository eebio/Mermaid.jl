using OrdinaryDiffEq, Random, Mermaid, BenchmarkTools, Plots, LinearAlgebra, Statistics
using Surrogates, Flux
using CellMLToolkit
# Random.seed!(1)

using Profile

ml = CellModel("examples/ohara_rudy_cipa_v1_2017.cellml")

tspan = (0, 1000.0)
longtime = 100000.0
n_samples = 1000
prob = ODEProblem(ml, tspan)
V = ml.sys.membrane₊v
Ki = ml.sys.intracellular_ions₊ki

# Slightly perturb initial conditions
prob[V] = -95.0
prob[Ki] = 144.0

comp1 = DEComponent(prob, Tsit5();
    name = "cell", time_step = 1000.0, state_names = Dict("v" => V, "ki" => Ki),
    intkwargs = (abstol = 1e-10, reltol = 1e-10, save_everystep = false, maxiters = Inf)
)

lb = prob.u0 - 0.01*abs.(prob.u0)
ub = prob.u0 + 0.01*abs.(prob.u0)

lb, ub = vec(lb - 0.1*(ub-lb)), vec(ub + 0.1*(ub-lb))

# Create a surrogate component
surrogate_comp = SurrogateComponent(comp1, lb, ub;
    n_samples = n_samples, n_epochs=10000)

println("=== Timing Comparison ===")

# Set up initial conditions for timing test
alg = MinimumTimeStepper()

# Benchmark the original component
println("Benchmarking original component...")
mp1 = MermaidProblem(components = [comp1], connectors = Connector[], max_t = longtime)
sol1 = solve(mp1, alg) # Compile first
ori_int = init(mp1, alg; save_vars = ["cell.v", "cell.ki"])
original_time = @benchmark solve!(int) setup = (int = deepcopy(ori_int))

# Benchmark the surrogate component
println("Benchmarking surrogate component...")
mp2 = MermaidProblem(components = [surrogate_comp], connectors = Connector[], max_t = longtime)
sol2 = solve(mp2, alg) # Compile first

surr_int = init(mp2, alg; save_vars = ["cell.v", "cell.ki"])
surrogate_time = @benchmark solve!(int) setup = (int = deepcopy(surr_int))

println("Original component minimum time: ", minimum(original_time.times) / 1e6, " ms")
println("Surrogate component minimum time: ", minimum(surrogate_time.times) / 1e6, " ms")
println("Speedup: ", minimum(original_time.times) / minimum(surrogate_time.times), "x")

println("\n=== Visualization ===")
# Plot comparison
# Plot the true time course over it
sol = solve(prob; tspan = longtime, abstol = 1e-10, reltol = 1e-10, maxiters = Inf)

p1 = plot(sol.t, sol[V], label = "True - Voltage", linewidth = 1)
p2 = plot(sol.t, sol[Ki], label = "True - Intracellular K+", linewidth = 1)

plot!(p1, sol1.t, sol1["cell.v"], label = "Original - Voltage", linewidth = 3)
plot!(p1, sol2.t, sol2["cell.v"], label = "Surrogate - Voltage", linestyle = :dash, linewidth = 3)
xlabel!(p1, "Time (ms)")
plot!(p2, sol1.t, sol1["cell.ki"], label = "Original - Intracellular K+", linewidth = 3)
plot!(p2, sol2.t, sol2["cell.ki"], label = "Surrogate - Intracellular K+", linestyle = :dash, linewidth = 3)
xlabel!(p2, "Time (ms)")
plot(p1, p2, layout = (2,1), size = (800,600))
ylims!(p1, -88.005, -87.92)
ylims!(p2, 144.475, 144.675)
ylabel!(p1, "Membrane Potential (mV)")
ylabel!(p2, "K+ Concentration (mM)")
plot(p1, p2, layout = (2, 1), size = (400, 500), top_margin = 3Plots.mm, right_margin=4Plots.mm, bottom_margin=-3Plots.mm, left_margin=1Plots.mm)
