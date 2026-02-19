"""
MIT License

Copyright (c) 2023-present The TrixiParticles.jl Authors (see AUTHORS.md)
Copyright (c) 2023-present Helmholtz-Zentrum hereon GmbH, Institute of Surface Science

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

@testitem "single trixi particles sim" begin
    # A portion of this code is adapted from the TrixiParticles.jl examples under the above
    # licence.

    # ==========================================================================================
    # 2D Hydrostatic Water Column Simulation
    #
    # This example simulates a column of water at rest in a tank under gravity.
    # It is a basic test case to verify hydrostatic pressure distribution and stability.
    # ==========================================================================================

    using TrixiParticles
    using OrdinaryDiffEq

    # ==========================================================================================
    # ==== Resolution
    fluid_particle_spacing = 0.05

    # Make sure that the kernel support of fluid particles at a boundary is always fully sampled
    boundary_layers = 3

    # ==========================================================================================
    # ==== Experiment Setup
    gravity = 9.81
    tspan = (0.0, 1.0)

    # Boundary geometry and initial fluid particle positions
    initial_fluid_size = (1.0, 0.9)
    tank_size = (1.0, 1.0)

    fluid_density = 1000.0
    sound_speed = 10.0
    state_equation = StateEquationCole(; sound_speed, reference_density = fluid_density,
        exponent = 7, clip_negative_pressure = false)

    tank = RectangularTank(
        fluid_particle_spacing, initial_fluid_size, tank_size, fluid_density,
        n_layers = boundary_layers, acceleration = (0.0, -gravity),
        state_equation = state_equation)

    # ==========================================================================================
    # ==== Fluid
    smoothing_length = 1.2 * fluid_particle_spacing
    smoothing_kernel = SchoenbergCubicSplineKernel{2}()

    alpha = 0.02
    viscosity_fluid = ArtificialViscosityMonaghan(alpha = alpha, beta = 0.0)

    fluid_density_calculator = ContinuityDensity()

    # This is to set acceleration with `trixi_include`
    system_acceleration = (0.0, -gravity)
    fluid_system = WeaklyCompressibleSPHSystem(tank.fluid, fluid_density_calculator,
        state_equation, smoothing_kernel,
        smoothing_length, viscosity = viscosity_fluid,
        acceleration = system_acceleration,
        source_terms = nothing)

    # ==========================================================================================
    # ==== Boundary

    # This is to set another boundary density calculation with `trixi_include`
    boundary_density_calculator = AdamiPressureExtrapolation()

    # This is to set wall viscosity with `trixi_include`
    viscosity_wall = nothing
    boundary_model = BoundaryModelDummyParticles(tank.boundary.density, tank.boundary.mass,
        state_equation = state_equation,
        boundary_density_calculator,
        smoothing_kernel, smoothing_length,
        viscosity = viscosity_wall)
    boundary_system = WallBoundarySystem(tank.boundary, boundary_model,
        prescribed_motion = nothing)

    # ==========================================================================================
    # ==== Simulation
    semi = Semidiscretization(fluid_system, boundary_system,
        parallelization_backend = PolyesterBackend())
    ode = semidiscretize(semi, tspan)

    # Use a Runge-Kutta method with automatic (error based) time step size control
    sol_trixi = solve(ode, RDPK3SpFSAL35(); saveat = 0.002, tstops = 0.002:0.002:1.0)

    comp = TrixiParticlesComponent(semi, RDPK3SpFSAL35(); name = "TrixiParticles Component", time_step = 0.002)
    mp = MermaidProblem(components = [comp], connectors = [], max_t = 1.0)
    alg = MinimumTimeStepper()
    sol_mermaid = solve(mp, alg; save_vars = ["TrixiParticles Component.#state"])
    a = [sol_mermaid(i)["TrixiParticles Component.#state"] for i in 0:0.1:1.0]
    b = [sol_trixi(i) for i in 0:0.1:1.0]
    @test a ≈ b
end

@testitem "test include" begin
    include(joinpath(@__DIR__, "included_file.jl"))
    @test a == 1
end

@testitem "coupled trixi particles sim" begin end

@testitem "state control" begin end

@testitem "interpolation" begin end
