using CommonSolve
using Agents
using SymbolicIndexingInterface

"""
    AgentsComponent <: AbstractTimeDependentComponent

A component that represents an agent-based model (ABM) using the Agents.jl package.

# Fields
- `model::StandardABM`: The agent-based model to be solved.
- `name::String="Agents Component"`: The name of the component.
- `state_names::Dict{String,Any} = Dict{String,Any}()`: A dictionary mapping ConnectedVariable names (as strings) to their corresponding properties in the agent model.
- `time_step::Float64=1.0`: The time step for the component (not the ABM solver timestep), i.e. how frequently should the inputs and outputs be updated.
"""
@kwdef struct AgentsComponent <: AbstractTimeDependentComponent
    model::StandardABM
    name::String = "Agents Component"
    state_names::Dict{String,Any} = Dict{String,Any}() # Naming model properties
    time_step::Float64 = 1.0 # TODO setup handling for time_step
end

"""
    AgentsComponentIntegrator <: ComponentIntegrator

A mutable struct that integrates an [AgentsComponent](@ref) using the Agents.jl package.

# Fields
- `integrator::StandardABM`: The agent-based model integrator.
- `component::AgentsComponent`: The AgentsComponent being integrated.
- `outputs::Dict{ConnectedVariable,Any}`: A dictionary mapping [ConnectedVariable](@ref) names to their initial values from the component.
- `inputs::Dict{ConnectedVariable,Any}`: A dictionary mapping [ConnectedVariable](@ref) names to their current values (initially 0).
"""
mutable struct AgentsComponentIntegrator <: ComponentIntegrator
    integrator::StandardABM
    component::AgentsComponent
    outputs::Dict{ConnectedVariable,Any}
    inputs::Dict{ConnectedVariable,Any}
end

"""
    init(c::AgentsComponent, conns::Vector{Connector})

Initializes an [AgentsComponentIntegrator](@ref) for the given [AgentsComponent](@ref) and its connections.

# Arguments
- `c::AgentsComponent`: The [AgentsComponent](@ref) to be integrated.
- `conns::Vector{Connector}`: The [Connectors](@ref Connector) that define the inputs and outputs of the component.

# Returns
- `AgentsComponentIntegrator`: The initialized integrator for the Agents.jl component.
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

Sets the state based on the current inputs and steps the [AgentsComponentIntegrator](@ref) in time one step.

# Arguments
- `compInt::AgentsComponentIntegrator`: The component integrator to be stepped. Its internal state will be mutated.
"""
function CommonSolve.step!(compInt::AgentsComponentIntegrator)
    for (key, value) in compInt.inputs
        setstate!(compInt, key, value)
    end
    step!(compInt.integrator)
end

"""
    getstate(compInt::AgentsComponentIntegrator, key::ConnectedVariable)

Retrieves the state of a specific variable from the [AgentsComponentIntegrator](@ref).

# Arguments
- `compInt::AgentsComponentIntegrator`: The component integrator containing the agent-based model.
- `key::ConnectedVariable`: The [ConnectedVariable](@ref) specifying which variable's state to retrieve.

# Returns
- The current state of the variable specified by `key`, which can be a model-level property or an agent-specific property.
"""
function getstate(compInt::AgentsComponentIntegrator, key::ConnectedVariable)
    if isnothing(key.variableindex)
        # No index for agent, so model-level property
        index = compInt.component.state_names[key.variable]
        return getproperty(compInt.integrator, index)
    elseif key.variableindex isa AbstractVector
        # Range of indexes for agents
        for i in key.variableindex
            # Index for agent, so get the agent property
            index = compInt.component.state_names[key.variable]
            return getproperty(compInt.integrator[i], index)
        end
    else
        # Single index for one agent
        index = compInt.component.state_names[key.variable]
        return getproperty(compInt.integrator[key.variableindex], index)
    end
end

"""
    setstate!(compInt::AgentsComponentIntegrator, key::ConnectedVariable, value)

Sets the state of a specific variable in the [AgentsComponentIntegrator](@ref).

# Arguments
- `compInt::AgentsComponentIntegrator`: The component integrator containing the agent-based model.
- `key::ConnectedVariable`: The [ConnectedVariable](@ref) specifying which variable's state to set.
- `value`: The value to assign to the specified variable's state.
"""
function setstate!(compInt::AgentsComponentIntegrator, key::ConnectedVariable, value)
    if isnothing(key.variableindex)
        # No index for agent, so model-level property
        index = compInt.component.state_names[key.variable]
        setproperty!(compInt.integrator, index, value)
    elseif key.variableindex isa AbstractVector
        # Range of indexes for agents
        for i in key.variableindex
            # Index for agent, so set the agent property
            index = compInt.component.state_names[key.variable]
            setproperty!(compInt.integrator[i], index, value)
        end
    else
        # Single index for one agent
        index = compInt.component.state_names[key.variable]
        setproperty!(compInt.integrator[key.variableindex], index, value)
    end
end

"""
    gettime(compInt::AgentsComponentIntegrator)

Returns the simulation time of the [AgentsComponentIntegrator](@ref) at the current state.

# Arguments
- `compInt::AgentsComponentIntegrator`: The component integrator for which to retrieve the current time.

# Returns
- `time`: The current simulation time.
"""
function gettime(compInt::AgentsComponentIntegrator)
    return abmtime(compInt.integrator)
end
