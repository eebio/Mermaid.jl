using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets, Plots, Mermaid
function get_nutrient_prob()
    # Parameters, variables, and derivatives
    @parameters t x y
    @parameters D γ
    @variables N(..) N₋(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    Dy = Differential(y)
    Dyy = Differential(y)^2

    # 1D PDE and boundary conditions
    eq = [
        Dt(N(t, x, y)) ~ 0.01*(Dxx(N(t, x, y)) + Dyy(N(t, x, y))) -
                    N₋(t, x, y),
        Dt(N₋(t, x, y)) ~ 0.0]
    bcs = [N(0, x, y) ~ x,
        N₋(0, x, y) ~ 0.0,
        Dx(N(t, 0, y)) ~ 0.0,
        Dx(N(t, 1, y)) ~ 0.0,
        Dx(N₋(t, 0, y)) ~ 0.0,
        Dx(N₋(t, 1, y)) ~ 0.0,
        Dy(N(t, x, 0)) ~ 0.0,
        Dy(N(t, x, 1)) ~ 0.0,
        Dy(N₋(t, x, 0)) ~ 0.0,
        Dy(N₋(t, x, 1)) ~ 0.0]

    # Space and time domains
    domains = [t ∈ Interval(0.0, 15.0),
        x ∈ Interval(0.0, 1.0),
        y ∈ Interval(0.0, 1.0)
    ]

    # PDE system
    @named pdesys = PDESystem(eq, bcs, domains, [t, x, y], [N(t, x, y), N₋(t, x, y)])

    # Method of lines discretization
    # Need a small dx here for accuracy
    dx = 0.1
    order = 2
    discretization = MOLFiniteDifference([x => dx, y => dx], t)

    # Convert the PDE problem into an ODE problem
    prob = discretize(pdesys, discretization);
    return prob
end
