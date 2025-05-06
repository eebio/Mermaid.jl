using DifferentialEquations, ModelingToolkit, MethodOfLines, DomainSets, Plots
# Parameters, variables, and derivatives
@parameters t x
@variables u(..) g(..) [irreducible=true]
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# 1D PDE and boundary conditions
eq = [Dt(u(t, x)) ~ Dxx(u(t, x)) + g(t, x),
    Dt(g(t, x)) ~ 1]
bcs = [u(0, x) ~ sin(pi * x),
    u(t, 0) ~ 0,
    u(t, 1) ~ 0,
    g(0, x) ~ 0.5,
    Dx(g(t, 0)) ~ 0,
    Dx(g(t, 1)) ~ 0]

# Space and time domains
domains = [t ∈ Interval(0.0, 1.0),
    x ∈ Interval(0.0, 1.0)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [g(t, x), u(t, x)])

dx = 0.1
# Method of lines discretization
discretization = MOLFiniteDifference([x => dx], t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)

solPDE = solve(prob, Euler(), dt=0.0001, adaptive=false)

display(plot(solPDE[u(t, x)]'))
display(plot(solPDE[g(t, x)]'))

# Can access u variables through this
ModelingToolkit.parse_variable(prob.f.sys, "u(t)")


# Define with Mermaid
using Mermaid

# Parameters, variables, and derivatives
# 1D PDE and boundary conditions
eq = [Dt(u(t, x)) ~ Dxx(u(t, x)) + g(t, x),
    Dt(g(t, x)) ~ 0]
bcs = [u(0, x) ~ sin(pi * x),
    u(t, 0) ~ 0,
    u(t, 1) ~ 0,
    g(0, x) ~ 0,
    Dx(g(t, 0)) ~ 0,
    Dx(g(t, 1)) ~ 0]

# Space and time domains
domains = [t ∈ Interval(0.0, 1.0),
    x ∈ Interval(0.0, 1.0)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [g(t, x), u(t, x)])

dx = 0.1
# Method of lines discretization
discretization = MOLFiniteDifference([x => dx], t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)

c1 = PDEComponent(
    model=prob,
    name="PDE",
    state_names=Dict(
        "u[2]" => 1, "u[3]" => 2, "u[4]" => 3, "u[5]" => 4, "u[6]" => 5,
        "u[7]" => 6, "u[8]" => 7, "u[9]" => 8, "u[10]" => 9,
        "g[2]" => 10, "g[3]" => 11, "g[4]" => 12,
     "g[5]" => 13, "g[6]" => 14, "g[7]" => 15, "g[8]" => 16, "g[9]" => 17,
      "g[10]" => 18), # Use hardcoded state indexes for now (12:22)
    time_step=0.0001,
    alg=Euler(),
    intkwargs=(:adaptive => false,),
)

function f2(u, p, t)
    return 1
end
u0 = 0.5
tspan = (0.0, 1.0)
prob = ODEProblem(f2, u0, tspan)
c2 = ODEComponent(
    model=prob,
    name="G",
    time_step=0.0001,
    state_names=Dict("g" => 1),
    alg=Euler(),
    intkwargs=(:adaptive => false,),
)

conn1 = Connector(
    inputs=["G.g"],
    outputs=["PDE.g[2]", "PDE.g[3]", "PDE.g[4]", "PDE.g[5]", "PDE.g[6]", "PDE.g[7]", "PDE.g[8]", "PDE.g[9]", "PDE.g[10]"],
)

conn2 = Connector(
    inputs=["PDE.u[2]", "PDE.u[3]", "PDE.u[4]", "PDE.u[5]", "PDE.u[6]", "PDE.u[7]", "PDE.u[8]", "PDE.u[9]", "PDE.u[10]"],
    outputs=[],
)

mp = MermaidProblem(components=[c1, c2], connectors=[conn1, conn2], max_t=1.0)

sol = solve(mp, MinimumTimeStepper())

# TODO putting incorrect names in connectors just skips then but should error
# TODO expand variable saving options: saveat, save_everystep, save_first, save_last, save_idx(s?),...
# TODO array connector (way of specifying a range of names easily)
finalsol = [0, sol.u["PDE.u[2]"][end], sol.u["PDE.u[3]"][end], sol.u["PDE.u[4]"][end], sol.u["PDE.u[5]"][end], sol.u["PDE.u[6]"][end], sol.u["PDE.u[7]"][end], sol.u["PDE.u[8]"][end], sol.u["PDE.u[9]"][end], sol.u["PDE.u[10]"][end], 0]
plot(finalsol, label="Mermaid")
plot!(solPDE[u(t,x)][end,:], label="MOL", linestyle=:dash)
