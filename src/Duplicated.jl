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
end

mutable struct DuplicatedComponentIntegrator <: ComponentIntegrator
    integrator::ComponentIntegrator
    component::DuplicatedComponent
    outputs::Dict{ConnectedVariable,Any}
    inputs::Dict{ConnectedVariable,Any}
    states::Vector
end

function CommonSolve.init(c::DuplicatedComponent, conns::Vector{Connector})
    integrator = init(c.component, conns)
    states = deepcopy(c.initial_states)
    inputs = deepcopy(integrator.inputs)
    outputs = deepcopy(integrator.outputs)
    # Values in inputs will be 0, but should be duplicated for each instance
    for key in keys(inputs)
        inputs[key] = zeros(c.instances)
    end
    # Values in outputs will be set based on the current state, but should instead match initial_states
    # Preallocate
    for key in keys(outputs)
        outputs[key] = zeros(c.instances)
    end
    for i in 1:c.instances
        # Set the state according to initial_states
        setstate!(integrator, states[i])
        for key in keys(outputs)
            newkey = ConnectedVariable(key.component, key.variable, nothing, nothing)
            outputs[key][i] = getstate(integrator, newkey)
        end
    end
    # Letting the internal integrator have inputs and outputs will break our setstate!
    integrator.inputs = Dict{ConnectedVariable, Any}()
    integrator.outputs = Dict{ConnectedVariable, Any}()
    # Create the DuplicatedComponentIntegrator
    integrator = DuplicatedComponentIntegrator(integrator, c, outputs, inputs, states)
    return integrator
end

function CommonSolve.step!(compInt::DuplicatedComponentIntegrator)
    t = gettime(compInt)
    for i in 1:compInt.component.instances
        # Set the time for this instance
        settime!(compInt.integrator, t)
        setstate!(compInt.integrator, compInt.states[i])
        # Set the inputs for this instance
        for (key, value) in compInt.inputs
            newkey = ConnectedVariable(key.component, key.variable, nothing, nothing)
            setstate!(compInt.integrator, newkey, value[i])
        end
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
    out = Vector{Any}(nothing, compInt.component.instances)
    if isnothing(key.duplicatedindex)
        index = 1:compInt.component.instances
    else
        index = key.duplicatedindex
    end
    for i in index
        # TODO I don't think this works, surely this doesnt change compInt.integrator?
        newkey = ConnectedVariable(key.component, key.variable, nothing, nothing)
        out[i] = getstate(compInt.integrator, newkey)
    end
    return out[index]
end

function setstate!(compInt::DuplicatedComponentIntegrator, key, value)
    for i in key.duplicatedindex
        setstate!(compInt.integrator, compInt.states[i])
        newkey = ConnectedVariable(key.component, key.variable, nothing, nothing)
        setstate!(compInt.integrator, newkey, value[i])
        states[i] = getstate(compInt.integrator, newkey)
    end
end

function gettime(compInt::DuplicatedComponentIntegrator)
    return gettime(compInt.integrator)
end
