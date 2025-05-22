using CommonSolve
using Agents
using SymbolicIndexingInterface

# Predefined concrete types
@kwdef struct AgentsComponent <: AbstractTimeDependentComponent
    model::StandardABM
    name::String = "Agents Component"
    state_names::Dict{String,Any} = Dict{String,Any}() # Naming model properties
    time_step::Float64 = 1.0 # TODO setup handling for time_step
end

mutable struct AgentsComponentIntegrator <: ComponentIntegrator
    integrator::StandardABM
    component::AgentsComponent
    outputs::Dict{ConnectedVariable,Any}
    inputs::Dict{ConnectedVariable,Any}
end

"""
    init(c::AgentsComponent)
"""
function CommonSolve.init(c::AgentsComponent, conns::Vector{Connector})
    outputs = Dict{ConnectedVariable,Any}() # Full variable name => Initial value from component
    inputs = Dict{ConnectedVariable,Any}() # Full variable name => Value (initially 0)
    for conn in conns
        # If connection has an input from this component, store its index and function as a ComponentIntegrator.output
        for input in conn.inputs
            if input.component == c.name
                outputProperty = c.state_names[split(input, ".")[2]]
                outputs[input] = getproperty(c.model, outputProperty)
            end
        end
        for output in conn.outputs
            if output.component == c.name
                inputs[output] = 0
            end
        end
    end

    integrator = AgentsComponentIntegrator(deepcopy(c.model), c, outputs, inputs)
    return integrator
end

"""
    step!(compInt::AgentsComponentIntegrator)

Steps the Agents.jl component integrator.
# Arguments
- `compInt::AgentsComponentIntegrator`: The Agents.jl component integrator to be stepped.
"""
function CommonSolve.step!(compInt::AgentsComponentIntegrator)
    for (key, value) in compInt.inputs
        setstate!(compInt, key, value)
    end
    step!(compInt.integrator)
end

function getstate(compInt::AgentsComponentIntegrator, key::ConnectedVariable)
    index = compInt.component.state_names[key.variable]
    return getproperty(compInt.integrator, index)
end

function setstate!(compInt::AgentsComponentIntegrator, key::ConnectedVariable, value)
    index = compInt.component.state_names[key.variable]
    setproperty!(compInt.integrator, index, value)
end

function gettime(compInt::AgentsComponentIntegrator)
    return abmtime(compInt.integrator)
end
