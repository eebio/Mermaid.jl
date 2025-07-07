using CommonSolve

"""
    DuplicatedComponent <: AbstractComponent

Represents a component that is duplicated in the simulation, allowing a single component to have multiple states.

# Fields
- `component::AbstractTimeDependentComponent`: The original component to be duplicated.
- `instances::Union{Int,Nothing}`: Number of instances of the component. If `nothing`, then the number is variable and determined by the simulation.
- `name::String`: Name of the duplicated component.
- `initial_states::Vector`: Vector of states for the duplicated component, where each state corresponds to a particular instance.
"""
@kwdef struct DuplicatedComponent <: AbstractTimeDependentComponent
    component::AbstractTimeDependentComponent
    instances::Union{Int,Nothing} = nothing
    name::String = component.name
    initial_states::Vector
    time_step::Float64 = component.time_step
    state_names = component.state_names
    default_state = zeros(length(initial_states[1]))
end

mutable struct DuplicatedComponentIntegrator <: ComponentIntegrator
    integrator::ComponentIntegrator
    component::DuplicatedComponent
    outputs::Dict{ConnectedVariable,Any}
    inputs::Dict{ConnectedVariable,Any}
    states::Vector
    ids::Union{Vector,Nothing}
end

function CommonSolve.init(c::DuplicatedComponent, conns::Vector{Connector})
    integrator = init(c.component, conns)
    states = deepcopy(c.initial_states)
    ids = isnothing(c.instances) ? [] : 1:c.instances

    outputs = Dict{ConnectedVariable,Any}() # Full variable name => Initial value from component
    inputs = Dict{ConnectedVariable,Any}() # Full variable name => Value (initially 0)
    for conn in conns
        # If connection has an input from this component, store its index and function as a ComponentIntegrator.output
        for input in conn.inputs
            if input.component == c.name
                outputs[input] = isnothing(c.instances) ? [] : zeros(c.instances)
            end
        end
        for output in conn.outputs
            if output.component == c.name
                inputs[output] = isnothing(c.instances) ? [] : zeros(c.instances)
            end
        end
    end

    # Values in outputs will be set based on the current state, but should instead match initial_states
    if !isnothing(c.instances)
        for i in 1:c.instances
            # Set the state according to initial_states
            setstate!(integrator, states[i])
            for key in keys(outputs)
                newkey = ConnectedVariable(key.component, key.variable, nothing, nothing)
                outputs[key][i] = getstate(integrator, newkey)
            end
        end
    end
    # Letting the internal integrator have inputs and outputs will break our setstate!
    integrator.inputs = Dict{ConnectedVariable, Any}()
    integrator.outputs = Dict{ConnectedVariable, Any}()
    # Create the DuplicatedComponentIntegrator
    integrator = DuplicatedComponentIntegrator(integrator, c, outputs, inputs, states, ids)
    return integrator
end

function CommonSolve.step!(compInt::DuplicatedComponentIntegrator)
    t = gettime(compInt)
    # Set the inputs for all states
    for (key, value) in compInt.inputs
        setstate!(compInt, key, value)
    end
    for i in eachindex(compInt.ids)
        # Set the time for this instance
        settime!(compInt.integrator, t)
        # Set the instances state
        setstate!(compInt.integrator, compInt.states[i])
        # Step the integrator for this instance
        step!(compInt.integrator)
        # Get the outputs for this instance
        for key in keys(compInt.outputs)
            newkey = ConnectedVariable(key.component, key.variable, nothing, nothing)
            compInt.outputs[key][i] = getstate(compInt.integrator, newkey)
        end
        # Store the state for this instance
        compInt.states[i] = getstate(compInt.integrator)
    end
end

function getstate(compInt::DuplicatedComponentIntegrator, key)
    if key.variable == "#ids"
        return compInt.ids
    end
    out = Vector{Any}(nothing, length(compInt.ids))
    if isnothing(key.duplicatedindex)
        index = eachindex(compInt.ids)
    else
        index = key.duplicatedindex
    end
    for i in index
        setstate!(compInt.integrator, compInt.states[i])
        newkey = ConnectedVariable(key.component, key.variable, key.variableindex, nothing)
        out[i] = getstate(compInt.integrator, newkey)
    end
    return out[index]
end

function setstate!(compInt::DuplicatedComponentIntegrator, key, value)
    if key.variable == "#ids"
        # Add any new states, remove any old states, and then set the ids vector
        states = similar(compInt.states, size(value))
        for i in eachindex(value)
            id = value[i]
            if id ∉ compInt.ids
                states[i] = deepcopy(compInt.component.default_state)
            else
                states[i] = compInt.states[findfirst(==(id), compInt.ids)]
            end
        end
        compInt.ids = value
        compInt.states = states
        for key in keys(compInt.outputs)
            resize!(compInt.outputs[key], length(compInt.ids))
        end
        for key in keys(compInt.inputs)
            resize!(compInt.inputs[key], length(compInt.ids))
        end
        return nothing
    end
    if isnothing(key.duplicatedindex)
        ids = compInt.ids
    else
        ids = key.duplicatedindex
    end
    for i in eachindex(ids)
        setstate!(compInt.integrator, compInt.states[i])
        newkey = ConnectedVariable(key.component, key.variable, nothing, nothing)
        setstate!(compInt.integrator, newkey, value[i])
        compInt.states[i] = getstate(compInt.integrator)
    end
end

function gettime(compInt::DuplicatedComponentIntegrator)
    return gettime(compInt.integrator)
end

function variables(component::DuplicatedComponent)
    return variables(component.component)
end
