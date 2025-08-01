@testitem "pde" begin
    using DifferentialEquations, ModelingToolkit, MethodOfLines, DomainSets
    # Parameters, variables, and derivatives
    @parameters t x
    @variables u(..) g(..) [irreducible = true]
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
        state_names=Dict("u" => [2:9..., 1], "g" => [11:18..., 10]), # TODO this needs to use the symbolic indexing
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
    c2 = DEComponent(
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
    finalsol = [0, sol(1)["PDE.u"]..., 0]
    @test all(isapprox.(finalsol, solPDE[u(t, x)][end, :], atol=1e-8))
end

@testitem "state control" begin
    using DifferentialEquations, ModelingToolkit, MethodOfLines, DomainSets
    @parameters t x
    @variables u(..) g(..) [irreducible = true]
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

    c1 = PDEComponent(
        model=prob,
        name="PDE",
        state_names=Dict("u" => [2:9..., 1], "g" => [11:18..., 10]), # TODO this needs to use the symbolic indexing
        time_step=0.01,
        alg=Tsit5(),
    )

    conn1 = Connector(
        inputs=["PDE.u[1:9]"],
        outputs=["other.u"],
    )
    conn2 = Connector(
        inputs=["PDE.g[1:9]"],
        outputs=["other.g"],
    )
    integrator = init(c1, [conn1, conn2])

    @test issetequal(Mermaid.variables(integrator), ["u", "g"])

    # Check initial state
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.u")) == [sin(pi * x) for x in 0.1:0.1:0.9]
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.g")) == [0.5 for _ in 0.1:0.1:0.9]
    @test Mermaid.getstate(integrator) == [[sin(pi * x) for x in [0.9, 0.1:0.1:0.8...]]...; [0.5 for _ in [0.9, 0.1:0.1:0.8...]]]
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.u[1]")) == sin(pi * 0.1)
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.g[1:3]")) == [0.5, 0.5, 0.5]
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.u[2:4]")) == [sin(pi * 0.2), sin(pi * 0.3), sin(pi * 0.4)]
    # Check setting state
    Mermaid.setstate!(integrator, ConnectedVariable("PDE.u"), [1.0 for _ in 0.1:0.1:0.9])
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.u")) == [1.0 for _ in 0.1:0.1:0.9]
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.g")) == [0.5 for _ in 0.1:0.1:0.9]
    Mermaid.setstate!(integrator, ConnectedVariable("PDE.u[1]"), 0.0)
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.u")) == [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    Mermaid.setstate!(integrator, ConnectedVariable("PDE.g[1:3]"), [0.0, 0.1, 0.2])
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.g")) == [0.0, 0.1, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]


    # Check time control
    @test Mermaid.gettime(integrator) == 0.0
    step!(integrator)
    @test Mermaid.gettime(integrator) == 0.01
    Mermaid.settime!(integrator, 0.1)
    @test Mermaid.gettime(integrator) == 0.1
    step!(integrator)
    @test Mermaid.gettime(integrator) == 0.11

    # Step means the state has changed
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.u")) ≠ [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.g")) ≠ [0.0, 0.1, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

    # Global setstate!
    Mermaid.setstate!(integrator, [-1.0 for _ in 1:18])
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.u")) == [-1.0 for _ in 0.1:0.1:0.9]
    @test Mermaid.getstate(integrator, ConnectedVariable("PDE.g")) == [-1.0 for _ in 0.1:0.1:0.9]

    # Test error on symbolic indexing
    c1 = PDEComponent(
        model=prob,
        name="PDE",
        state_names=Dict("u" => u, "g" => g),
        time_step=0.01,
        alg=Tsit5(),
    )
    @test_throws ArgumentError init(c1, [conn1, conn2])
end
