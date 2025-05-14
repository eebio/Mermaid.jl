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
    outputs::Dict{String,Any}
    inputs::Dict{String,Any}
end

"""
    init(c::AgentsComponent)
"""
function CommonSolve.init(c::AgentsComponent, conns::Vector{Connector})
    outputs = Dict{String,Any}() # Full variable name => Initial value from component
    inputs = Dict{String,Any}() # Full variable name => Value (initially 0)
    for conn in conns
        # If connection has an input from this component, store its index and function as a ComponentIntegrator.output
        for input in conn.inputs
            if split(input, ".")[1] != c.name
                continue
            end
            outputProperty = c.state_names[split(input, ".")[2]]

            # TODO change for agents
            outputs[input] = getproperty(c.model, outputProperty)
        end
        for output in conn.outputs
            if split(output, ".")[1] != c.name
                continue
            end
            inputs[output] = 0
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
        index = compInt.component.state_names[split(key, ".")[2]]
        setproperty(compInt.integrator.model, index, value)
    end
    step!(compInt.integrator)
end

function getstate(compInt::AgentsComponentIntegrator, state_name::String)
    index = compInt.component.state_names[state_name]
    return getproperty(compInt.integrator.model, index)
end

function setstate!(compInt::AgentsComponentIntegrator, state_name::String, value)
    index = compInt.component.state_names[state_name]
    setproperty(compInt.integrator.model, index, value)
end

function gettime(compInt::AgentsComponentIntegrator)
    return abmtime(compInt.integrator)
end
