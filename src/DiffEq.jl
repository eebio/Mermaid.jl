using CommonSolve

"""
    init(c::ODEComponent)
"""
function CommonSolve.init(c::ODEComponent)
    # Initialize the integrator for the component
    outputs = Dict{String,Any}([key => c.model.u0[value] for (key, value) in c.output_indices])
    inputs = Dict{String,Any}([key => 0.0 for key in c.input_names])
    integrator = ODEComponentIntegrator(init(c.model, c.alg; dt=c.time_step, c.intkwargs...), c, outputs, inputs)
    return integrator
end

"""
    step!(compInt::ODEComponentIntegrator)

Steps the ODE component integrator.
# Arguments
- `compInt::ODEComponentIntegrator`: The ODE component integrator to be stepped.
"""
function CommonSolve.step!(compInt::ODEComponentIntegrator)
    # Update compInt.integrator with the parameters from inputs
    p = collect(values(compInt.inputs)) # Need to be smarter about this
    compInt.integrator.p = p
    u_modified!(compInt.integrator, true)
    CommonSolve.step!(compInt.integrator)
end
