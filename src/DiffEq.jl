using ModelingToolkit
using DifferentialEquations
using OrdinaryDiffEqCore
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

@kwdef mutable struct ODEComponentIntegrator <: ComponentIntegrator
    integrator::OrdinaryDiffEqCore.ODEIntegrator
    component::ODEComponent
    outputs::Dict{ConnectedVariable,Any} = Dict{ConnectedVariable,Any}()
    inputs::Dict{ConnectedVariable,Any} = Dict{ConnectedVariable,Any}()
end

function CommonSolve.init(c::ODEComponent, conns::Vector{Connector})
    integrator = ODEComponentIntegrator(integrator = init(c.model, c.alg; dt=c.time_step, c.intkwargs...), component = c)
    inputs, outputs = inputsandoutputs(integrator, conns)
    integrator.inputs = inputs
    integrator.outputs = outputs
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

function variables(component::ODEComponent)
    return keys(component.state_names)
end
