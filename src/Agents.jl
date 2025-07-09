using Agents

"""
    AgentsComponent <: AbstractTimeDependentComponent

A component that represents an agent-based model (ABM) using the Agents.jl package.

# Fields
- `model::StandardABM`: The agent-based model to be solved.
- `name::String="Agents Component"`: The name of the component.
- `state_names::Dict{String,Any} = Dict{String,Any}()`: A dictionary mapping ConnectedVariable names (as strings) to their corresponding properties in the agent model. The properties differentiate between `:model` properties and `:agent` properties.
- `time_step::Float64=1.0`: The time step for the component (not the ABM solver timestep), i.e. how frequently should the inputs and outputs be updated.
"""
@kwdef struct AgentsComponent <: AbstractTimeDependentComponent
    model::StandardABM
    name::String = "Agents Component"
    state_names::Dict{String,Any} = Dict{String,Any}() # Naming model properties
    time_step::Float64 = 1.0
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
@kwdef mutable struct AgentsComponentIntegrator <: ComponentIntegrator
    integrator::StandardABM
    component::AgentsComponent
    outputs::OrderedDict{ConnectedVariable,Any} = OrderedDict{ConnectedVariable,Any}()
    inputs::OrderedDict{ConnectedVariable,Any} = OrderedDict{ConnectedVariable,Any}()
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
    integrator = AgentsComponentIntegrator(integrator = deepcopy(c.model), component = c)
    inputs, outputs = inputsandoutputs(integrator, conns)
    integrator.inputs = inputs
    integrator.outputs = outputs
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
    if key.variable == "#model"
        return getstate(compInt, true)
    elseif key.variable == "#ids"
        return collect(allids(compInt.integrator))
    end
    index = compInt.component.state_names[key.variable]
    if isnothing(key.variableindex)
        # If model level property exists, return it directly
        !isnothing(abmproperties(compInt.integrator)) && haskey(abmproperties(compInt.integrator), index) && return getproperty(compInt.integrator, index)
        # Otherwise, assume it's an agent property and return it for all agents
        return [getproperty(i, index) for i in allagents(compInt.integrator)]
    else
        # If model level property exists, return it after indexing
        !isnothing(abmproperties(compInt.integrator)) && haskey(abmproperties(compInt.integrator), index) && return getproperty(compInt.integrator, index)[key.variableindex]
        # Otherwise, assume it's an agent property and return it for all agents in the specified range
        out = [getproperty(compInt.integrator[i], index) for i in key.variableindex]
        if length(out) == 1
            return out[1] # If only one agent, return the single value
        else
            return out # Otherwise, return the vector of values
        end
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
    # TODO add ids exception? Is there a way to specify the id when creating an agent
    if key.variable == "#model"
        setstate!(compInt, value)
        return nothing
    end
    index = compInt.component.state_names[key.variable]
    if isnothing(key.variableindex)
        # If model level property exists, return it directly
        if !isnothing(abmproperties(compInt.integrator)) && haskey(abmproperties(compInt.integrator), index)
            setindex!(abmproperties(compInt.integrator), value, index)
        else
            # Otherwise, assume it's an agent property and set it for all agents
            for (i, agent) in enumerate(allagents(compInt.integrator))
                setproperty!(agent, index, value[i])
            end
        end
    else
        # If model level property exists, return it after indexing
        if !isnothing(abmproperties(compInt.integrator)) && haskey(abmproperties(compInt.integrator), index)
            k = 1
            for i in key.variableindex
                abmproperties(compInt.integrator)[index][i] = value[k]
                k += 1
            end
        else
            # Otherwise, assume it's an agent property and set it for all agents in the specified range
            k = 1
            for i in key.variableindex
                setproperty!(compInt.integrator[i], index, value[k])
                k += 1
            end
        end
    end
    return nothing
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
    return abmtime(compInt.integrator)*compInt.component.time_step
end

"""
    getstate(compInt::AgentsComponentIntegrator)

Returns the current state of the [AgentsComponentIntegrator](@ref).

# Arguments
- `compInt::AgentsComponentIntegrator`: The component integrator for which to retrieve the current state.
- `copy::Bool=false`: If true, returns a deep copy of the state.

# Returns
- `state::StandardABM`: The current state of the agent-based model being integrated.
"""
function getstate(compInt::AgentsComponentIntegrator, copy::Bool=false)
    if copy
        return deepcopy(compInt.integrator)
    else
        return compInt.integrator
    end
end

"""
    setstate!(compInt::AgentsComponentIntegrator, state::StandardABM)

Sets the state of the [AgentsComponentIntegrator](@ref) to a new state.

# Arguments
- `compInt::AgentsComponentIntegrator`: The component integrator to be updated.
- `state::StandardABM`: The new state to set for the component integrator.
"""
function setstate!(compInt::AgentsComponentIntegrator, state::StandardABM)
    compInt.integrator = state
end

function settime!(compInt::AgentsComponentIntegrator, time::Float64)
    # Agents.jl does not support setting time directly, but the time is stored within the state of the integrator.
    return nothing
end

function variables(component::AgentsComponent)
    return union(keys(component.state_names), ["#model"])
end
