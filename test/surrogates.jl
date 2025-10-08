@testitem "surrogate integrator" begin
    using DifferentialEquations
    # Define a simple ODE: dx/dt = -x, x(0) = 1
    f(u, p, t) = -u
    u0 = 1.0
    tspan = (0.0, 1.0)
    prob = ODEProblem(f, u0, tspan)
    state_names = Dict("x" => 1)
    ode_comp = DEComponent(model=prob, name="decay", state_names=state_names, time_step=0.1)

    # Set bounds for surrogate sampling
    lower = [0.0]
    upper = [1.0]

    # Create surrogate component
    surrogate_comp = SurrogateComponent(
        component=ode_comp,
        lower_bound=lower,
        upper_bound=upper,
    )

    # Dummy connector list (no connections for this test)
    conns = Connector[]

    # Initialize both components
    ode_int = init(ode_comp, conns)
    surrogate_int = init(surrogate_comp, conns)

    # Step both integrators and compare
    n_steps = 20
    tol = 0.01  # Surrogate tolerance

    for i in 1:n_steps
        step!(ode_int)
        step!(surrogate_int)
        orig = Mermaid.getstate(ode_int, ConnectedVariable("decay.x"))
        sur = Mermaid.getstate(surrogate_int, ConnectedVariable("decay.x"))
        @test isapprox(sur, orig; atol=tol)
    end
end

@testitem "surrogate integrator" begin
    using DifferentialEquations
    # Define a simple ODE: dx/dt = -x, x(0) = 1
    f(u, p, t) = [-u[1], u[2]-u[1]]
    u0 = [1.0, 0.5]
    tspan = (0.0, 1.0)
    prob = ODEProblem(f, u0, tspan)
    state_names = Dict("x" => 1)
    ode_comp = DEComponent(model=prob, name="decay", state_names=state_names, time_step=0.1)

    # Set bounds for surrogate sampling
    lower = [0.0, 0.0]
    upper = [1.0, 1.0]

    # Create surrogate component
    surrogate_comp = SurrogateComponent(
        component=ode_comp,
        lower_bound=lower,
        upper_bound=upper,
        n_samples=2000,
        n_epochs=2000,
    )

    # Dummy connector list (no connections for this test)
    conns = Connector[]

    # Initialize both components
    ode_int = init(ode_comp, conns)
    surrogate_int = init(surrogate_comp, conns)

    # Step both integrators and compare
    n_steps = 20
    tol = 0.1  # Looser tolerances to save training time in tests

    for i in 1:n_steps
        step!(ode_int)
        step!(surrogate_int)
        orig = Mermaid.getstate(ode_int, ConnectedVariable("decay.x"))
        sur = Mermaid.getstate(surrogate_int, ConnectedVariable("decay.x"))
        @test isapprox(sur, orig; atol=tol)
    end
end
