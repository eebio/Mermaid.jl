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

using SymbolicIndexingInterface
function var_index(s)
    return s => variable_index(prob, ModelingToolkit.parse_variable(prob.f.sys, s))
end

c1 = PDEComponent(
    model=prob,
    name="PDE",
    state_names=Dict("u" => [2:9..., 1], "g" => [11:18..., 10]),
    time_step=0.0001,
    alg=Euler(),
    intkwargs=(:adaptive => false,),
)

function f2(u, p, t)
    return 1.0
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

conn = Connector(
    inputs=["G.g"],
    outputs=["PDE.g[1:9]"],
)

mp = MermaidProblem(components=[c1, c2], connectors=[conn], max_t=1.0)
sol = solve(mp, MinimumTimeStepper())
@show sol["PDE.u"][end-1]
@show solPDE[u(t, x)][end, :]

# TODO expand variable saving options: saveat, save_everystep, save_first, save_last, save_idx(s?),...
finalsol = [0, sol["PDE.u"][end-1]..., 0]
plot(finalsol, label="Mermaid")
plot!(solPDE[u(t,x)][end,:], label="MOL", linestyle=:dash)
