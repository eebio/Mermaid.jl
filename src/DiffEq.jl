using CommonSolve
using ModelingToolkit
using SymbolicIndexingInterface

"""
    ODEComponent <: AbstractTimeDependentComponent
A component that represents an ODE system, defined by an [ODEProblem](@extref DiffEq types/ode_types).

# Fields
- `model::ODEProblem`: The ODE problem to be solved.
- `name::String="ODE Component"`: The name of the component.
- `state_names::Dict{String,Any} = Dict{String,Any}()`: A dictionary mapping variable names (as strings) to their corresponding indices in the state vector or symbols from [ModelingToolkit](@extref ModelingToolkit index)/[Symbolics](@extref Symbolics index).
- `time_step::Float64=1.0`: The time step for the component (not the ODE solver timestep), i.e. how frequently should the inputs and outputs be updated.
- `alg=Rodas5()`: The algorithm used for solving the ODEProblem.
- `intkwargs::Tuple=()`: Additional keyword arguments for the ODE solver.
"""
@kwdef struct ODEComponent <: AbstractTimeDependentComponent
    model::ODEProblem
    name::String = "ODE Component"
    state_names::Dict{String,Any} = Dict{String,Any}()
    time_step::Float64 = 1.0
    alg = Rodas5()
    intkwargs::Tuple = ()
end

mutable struct ODEComponentIntegrator <: ComponentIntegrator
    integrator::OrdinaryDiffEqCore.ODEIntegrator
    component::ODEComponent
    outputs::Dict{ConnectedVariable,Any}
    inputs::Dict{ConnectedVariable,Any}
end

function CommonSolve.init(c::ODEComponent, conns::Vector{Connector})
    outputs = Dict{ConnectedVariable, Any}() # Full variable name => Initial value from component
    inputs = Dict{ConnectedVariable, Any}() # Full variable name => Value (initially 0)
    for conn in conns
        # If connection has an input from this component, store its index and function as a ComponentIntegrator.output
        for input in conn.inputs
            if input.component == c.name
                outputIndex = c.state_names[input.variable]
                # If the index is a MTK symbol then get the variable index
                if symbolic_type(outputIndex) != NotSymbolic()
                    outputIndex = variable_index(c.model.f.sys, outputIndex)
                end
                # TODO I think we can allow variable indexes here, for if each element of the state is a vector
                outputs[input] = c.model.u0[outputIndex]
            end
        end
        for output in conn.outputs
            if output.component == c.name
                inputs[output] = 0
            end
        end
    end

    integrator = ODEComponentIntegrator(init(c.model, c.alg; dt=c.time_step, c.intkwargs...), c, outputs, inputs)
    return integrator
end

function CommonSolve.step!(compInt::ODEComponentIntegrator)
    for (key, value) in compInt.inputs
        setstate!(compInt, key, value)
    end
    u_modified!(compInt.integrator, true)
    CommonSolve.step!(compInt.integrator, compInt.component.time_step, true)
end

function getstate(compInt::ODEComponentIntegrator, key)
    index = compInt.component.state_names[key.variable]
    return compInt.integrator[index]
end

function getstate(compInt::ODEComponentIntegrator)
    return compInt.integrator.u
end

function setstate!(compInt::ODEComponentIntegrator, key, value)
    index = compInt.component.state_names[key.variable]
    compInt.integrator[index] = value
end

function setstate!(compInt::ODEComponentIntegrator, value)
    compInt.integrator.u = value
end

function gettime(compInt::ODEComponentIntegrator)
    return compInt.integrator.t
end

function settime!(compInt::ODEComponentIntegrator, t)
    u_modified!(compInt.integrator, true)
    compInt.integrator.t = t
end
