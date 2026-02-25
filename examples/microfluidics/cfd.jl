using TrixiParticles
using OrdinaryDiffEq
using Plots

tspan = (0.0, 2.0)

block_width = 1.0
height = 1.0
reynolds_number = 100
fluid_density = 1000.0
flow_direction = [1.0, 0.0]

const prescribed_velocity = (1.0, 0.0)

sound_speed = 10 * maximum(abs.(prescribed_velocity))

particle_spacing = 0.05
boundary_layers = 4
open_boundary_layers = 6
boundary_thickness = boundary_layers * particle_spacing
open_boundary_size = (particle_spacing * open_boundary_layers, height)

# Create from a bunch of tanks
tank1 = RectangularTank(particle_spacing, (block_width, height),
    (block_width, height), fluid_density,
    n_layers = boundary_layers, velocity = prescribed_velocity,
    faces = (false, false, true, true),
    coordinates_eltype = Float64)

tank2 = RectangularTank(particle_spacing, (block_width, height),
    (block_width, height), fluid_density,
    n_layers = boundary_layers,
    min_coordinates = (block_width, -height),
    velocity = prescribed_velocity,
    faces = (true, true, true, false),
    coordinates_eltype = Float64)

tank3 = RectangularTank(particle_spacing, (block_width, height),
    (block_width, height), fluid_density,
    n_layers = boundary_layers,
    min_coordinates = (block_width, 0.0),
    velocity = prescribed_velocity,
    faces = (false, false, false, true),
    coordinates_eltype = Float64)

tank4 = RectangularTank(particle_spacing, (block_width, height),
    (block_width, height), fluid_density,
    n_layers = boundary_layers,
    min_coordinates = (2*block_width, 0.0),
    velocity = prescribed_velocity,
    faces = (false, false, true, true),
    coordinates_eltype = Float64)

min_coords_inlet = (-open_boundary_layers * particle_spacing, 0.0)
inlet = RectangularTank(particle_spacing, open_boundary_size, open_boundary_size,
    fluid_density, n_layers = boundary_layers,
    min_coordinates = min_coords_inlet,
    faces = (false, false, true, true),
    coordinates_eltype = Float64)

min_coords_outlet = (3*block_width, 0.0)
outlet = RectangularTank(particle_spacing, open_boundary_size, open_boundary_size,
    fluid_density, n_layers = boundary_layers,
    min_coordinates = min_coords_outlet,
    faces = (false, false, true, true),
    coordinates_eltype = Float64)

fluid = union(tank1.fluid, tank2.fluid, tank3.fluid, tank4.fluid)
boundary = union(tank1.boundary, tank2.boundary, tank3.boundary, tank4.boundary)
plot(fluid, boundary)
plot!(inlet.fluid, inlet.boundary)
plot!(outlet.fluid, outlet.boundary)

# ==========================================================================================
# ==== Fluid

smoothing_length = 1.5 * particle_spacing
smoothing_kernel = WendlandC2Kernel{2}()

fluid_density_calculator = ContinuityDensity()

kinematic_viscosity = maximum(prescribed_velocity) * height / reynolds_number
n_buffer_particles = convert(Int, 10 * height/particle_spacing)
viscosity = ViscosityAdami(nu = kinematic_viscosity)

# Alternatively the WCSPH scheme can be used
state_equation = StateEquationCole(; sound_speed, reference_density = fluid_density,
    exponent = 1)
density_diffusion = DensityDiffusionMolteniColagrossi(delta = 0.1)

fluid_system = WeaklyCompressibleSPHSystem(fluid, fluid_density_calculator,
    state_equation, smoothing_kernel,
    density_diffusion = density_diffusion,
    smoothing_length, viscosity = viscosity,
    shifting_technique = ParticleShiftingTechnique(v_max_factor = 1.5),
    buffer_size = n_buffer_particles)

function velocity_function2d(pos, t)
    # Use this for a time-dependent inflow velocity
    # return SVector(0.5prescribed_velocity * sin(2pi * t) + prescribed_velocity, 0)

    return SVector(prescribed_velocity)
end

open_boundary_model = BoundaryModelMirroringTafuni(; mirror_method = ZerothOrderMirroring())

reference_velocity_in = velocity_function2d
reference_pressure_in = nothing
reference_density_in = nothing
boundary_type_in = InFlow()
face_in = ([0.0, 0.0], [0.0, height])
inflow = BoundaryZone(; boundary_face = face_in, face_normal = flow_direction,
    open_boundary_layers, density = fluid_density, particle_spacing,
    reference_density = reference_density_in,
    reference_pressure = reference_pressure_in,
    reference_velocity = reference_velocity_in,
    initial_condition = inlet.fluid, boundary_type = boundary_type_in)

reference_velocity_out = nothing
reference_pressure_out = nothing
reference_density_out = nothing
boundary_type_out = OutFlow()
face_out = ([3*block_width, 0.0], [3*block_width, height])
outflow = BoundaryZone(; boundary_face = face_out, face_normal = (-flow_direction),
    open_boundary_layers, density = fluid_density, particle_spacing,
    reference_density = reference_density_out,
    reference_pressure = reference_pressure_out,
    reference_velocity = reference_velocity_out,
    initial_condition = outlet.fluid, boundary_type = boundary_type_out)

open_boundary = OpenBoundarySystem(inflow, outflow; fluid_system,
    boundary_model = open_boundary_model,
    buffer_size = n_buffer_particles)

# ==========================================================================================
# ==== Boundary
wall = union(boundary, inlet.boundary, outlet.boundary)
viscosity_boundary = viscosity
boundary_model = BoundaryModelDummyParticles(wall.density, wall.mass,
    AdamiPressureExtrapolation(),
    state_equation = state_equation,
    viscosity = viscosity_boundary,
    smoothing_kernel, smoothing_length)

boundary_system = WallBoundarySystem(wall, boundary_model)

# ==========================================================================================
# ==== Simulation
min_corner = minimum(wall.coordinates .- particle_spacing, dims = 2)
max_corner = maximum(wall.coordinates .+ particle_spacing, dims = 2)

nhs = GridNeighborhoodSearch{2}(;
    cell_list = FullGridCellList(; min_corner, max_corner),
    update_strategy = ParallelUpdate())

semi = Semidiscretization(fluid_system, open_boundary, boundary_system,
    neighborhood_search = nhs,
    parallelization_backend = PolyesterBackend())

ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval = 100)
saving_callback = SolutionSavingCallback(;dt = 0.02, prefix = "")

extra_callback = nothing

callbacks = CallbackSet(info_callback, saving_callback, UpdateCallback(), extra_callback)

sol = solve(ode, RDPK3SpFSAL35(),
    abstol = 1e-5, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
    reltol = 1e-3, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
    dtmax = 1e-2, # Limit stepsize to prevent crashing
    save_everystep = false, callback = callbacks);
