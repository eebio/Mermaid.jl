using CommonSolve
using ModelingToolkit
using SymbolicIndexingInterface

"""
    init(c::ODEComponent)
"""
function CommonSolve.init(c::ODEComponent)
    # Initialize the integrator for the component
    outputs = Dict{String,Any}([key => c.model.u0[value] for (key, value) in c.output_indices])
    inputs = Dict{String,Any}([key => 0.0 for key in c.input_names])

    if isa(c.model.p, MTKParameters)
        setinputs! = setp(c.model.f.sys, [c.mtk_input_symbols[names] for names in c.input_names])
    else
        function setinputs!(int::ODEIntegrator, inputs)
            int.p[1:length(inputs)] = inputs
        end
    end

    integrator = ODEComponentIntegrator(init(c.model, c.alg; dt=c.time_step, c.intkwargs...), c, outputs, inputs, setinputs!)
    return integrator
end

"""
    step!(compInt::ODEComponentIntegrator)

Steps the ODE component integrator.
# Arguments
- `compInt::ODEComponentIntegrator`: The ODE component integrator to be stepped.
"""
function CommonSolve.step!(compInt::ODEComponentIntegrator)
    compInt.setinputs!(compInt.integrator, collect(values(compInt.inputs)))
    u_modified!(compInt.integrator, true)
    CommonSolve.step!(compInt.integrator)
end
