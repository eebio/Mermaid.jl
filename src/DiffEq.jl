using ModelingToolkit
using SciMLBase
using DifferentialEquations
using OrdinaryDiffEqCore
using SymbolicIndexingInterface

"""
    DEComponent <: AbstractTimeDependentComponent
A component that represents an DE system. For example, defined by an [ODEProblem](@extref DiffEq types/ode_types).

# Fields
- `model::SciML.AbstractDEProblem`: The SciML DE problem to be solved.
- `name::String="DE Component"`: The name of the component.
- `state_names::Dict{String,Any} = Dict{String,Any}()`: A dictionary mapping variable names (as strings) to their corresponding indices in the state vector or symbols from [ModelingToolkit](@extref ModelingToolkit index)/[Symbolics](@extref Symbolics index).
- `time_step::Float64=1.0`: The time step for the component (not the DE solver timestep), i.e. how frequently should the inputs and outputs be updated.
- `alg=Rodas5()`: The algorithm used for solving the DEProblem.
- `intkwargs::Tuple=()`: Additional keyword arguments for the DE solver.
"""
@kwdef struct DEComponent <: AbstractTimeDependentComponent
    model::SciMLBase.AbstractDEProblem
    name::String = "DE Component"
    state_names::Dict{String,Any} = Dict{String,Any}()
    time_step::Float64 = 1.0
    alg = Rodas5()
    intkwargs::Tuple = ()
end

@kwdef mutable struct DEComponentIntegrator <: ComponentIntegrator
    integrator::SciMLBase.DEIntegrator
    component::DEComponent
    outputs::OrderedDict{ConnectedVariable,Any} = OrderedDict{ConnectedVariable,Any}()
    inputs::OrderedDict{ConnectedVariable,Any} = OrderedDict{ConnectedVariable,Any}()
end

function CommonSolve.init(c::DEComponent, conns::Vector{Connector})
    integrator = DEComponentIntegrator(integrator = init(c.model, c.alg; dt=c.time_step, c.intkwargs...), component = c)
    inputs, outputs = inputsandoutputs(integrator, conns)
    integrator.inputs = inputs
    integrator.outputs = outputs
    return integrator
end

function CommonSolve.step!(compInt::DEComponentIntegrator)
    for (key, value) in compInt.inputs
        setstate!(compInt, key, value)
    end
    CommonSolve.step!(compInt.integrator, compInt.component.time_step, true)
end

function getstate(compInt::DEComponentIntegrator, key)
    if first(key.variable) == '#'
        if key.variable == "#time"
            return compInt.integrator.t
        end
    end
    index = compInt.component.state_names[key.variable]
    return compInt.integrator[index]
end

function getstate(compInt::DEComponentIntegrator)
    return compInt.integrator.u
end

function setstate!(compInt::DEComponentIntegrator, key, value)
    u_modified!(compInt.integrator, true)
    if first(key.variable) == '#'
        if key.variable == "#time"
            compInt.integrator.t = value
            return nothing
        end
    end
    index = compInt.component.state_names[key.variable]
    compInt.integrator[index] = value
end

function setstate!(compInt::DEComponentIntegrator, value)
    u_modified!(compInt.integrator, true)
    compInt.integrator.u = value
end

function variables(component::DEComponent)
    return keys(component.state_names)
end
