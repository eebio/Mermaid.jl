"""
    DuplicatedComponent <: AbstractComponent

Represents a component that is duplicated in the simulation, allowing a single component to
    have multiple states.

# Fields
- `component::AbstractTimeDependentComponent`: The original component to be duplicated.
- `instances::Union{Int,Nothing}`: Number of instances of the component. If `nothing`, then
    the number is variable and determined by the simulation.
- `name::String`: Name of the duplicated component.
- `init_states::Vector`: Vector of states for the duplicated component, where each state
    corresponds to a particular instance.
"""
struct DuplicatedComponent{T, U, V, W} <:
       AbstractTimeDependentComponent where {T <: AbstractTimeDependentComponent}
    component::T
    instances::Union{Int, Nothing}
    name::String
    init_states::U
    time_step::Float64
    state_names::V
    default_state::W
end

"""
    DuplicatedComponent(args...; kwargs...)

Duplicate an existing component into a [DuplicatedComponent](@ref).

# Arguments
- `component::AbstractTimeDependentComponent`: The component to be duplicated.
- `init_states::AbstractVector`: A vector of initial states for the duplicated component.

# Keyword Arguments
- `instances::Union{Int, Nothing}`: Number of instances of the component. If `nothing`, then
    the number is variable and determined by the simulation. Defaults to `nothing`.
- `name::AbstractString`: Name of the duplicated component. Defaults to the original
    component's name.
- `time_step::Real`: Time step for the duplicated component. Defaults to the original
    component's time step.
- `state_names`: A dictionary mapping variable names (as strings) to their corresponding
    variables in the original component. Defaults to the original component's state names.
- `default_state`: The default state to use when a new instance is created. Defaults to a
    zero vector of the same length as the first initial state.
"""
function DuplicatedComponent(component::AbstractTimeDependentComponent,
        init_states::AbstractVector; instances::Union{Int, Nothing} = nothing,
        name::AbstractString = component.name, time_step::Real = component.time_step,
        state_names = component.state_names,
        default_state = zeros(length(init_states[1])))
    return DuplicatedComponent(component, instances, name, init_states, time_step,
        state_names, default_state)
end

mutable struct DuplicatedComponentIntegrator{T <: AbstractComponentIntegrator} <:
               AbstractComponentIntegrator
    integrator::T
    component::DuplicatedComponent
    states::Vector
    ids::Vector{Int}
    init_states::Dict
end

function CommonSolve.init(c::DuplicatedComponent)
    integrator = CommonSolve.init(c.component)
    states = deepcopy(c.init_states)
    ids = isnothing(c.instances) ? [] : 1:(c.instances)

    # Create the DuplicatedComponentIntegrator
    return DuplicatedComponentIntegrator{typeof(integrator)}(
        integrator, c, states, ids, Dict())
end

function CommonSolve.step!(compInt::DuplicatedComponentIntegrator)
    t = gettime(compInt)
    for i in eachindex(compInt.ids)
        # Set the time for this instance
        settime!(compInt, t)
        # Set the instances state
        setstate!(compInt.integrator, compInt.states[i])
        # Step the integrator for this instance
        CommonSolve.step!(compInt.integrator)
        # Store the state for this instance
        compInt.states[i] = getstate(compInt.integrator; copy = true)
    end
end

function getstate(compInt::DuplicatedComponentIntegrator, key)
    if first(key.variable) == '#'
        if key.variable == "#time"
            return getstate(compInt.integrator, key)
        end
        if key.variable == "#ids"
            return compInt.ids
        end
        if key.variable == "#states"
            return compInt.states
        end
        if key.variable == "#init_states"
            return compInt.init_states
        end
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
    if first(key.variable) == '#'
        if key.variable == "#time"
            setstate!(compInt.integrator, key, value)
            return nothing
        end
        if key.variable == "#ids"
            # Add any new states, remove any old states, and then set the ids vector
            states = similar(compInt.states, size(value))
            for i in eachindex(value)
                id = value[i]
                if id ∉ compInt.ids
                    if id ∈ keys(compInt.init_states)
                        states[i] = deepcopy(compInt.init_states[id])
                    else
                        # If no initial state is provided, use the default state
                        states[i] = deepcopy(compInt.component.default_state)
                    end
                else
                    states[i] = compInt.states[findfirst(==(id), compInt.ids)]
                end
            end
            compInt.ids = value
            compInt.states = states
            return nothing
        end
        if key.variable == "#states"
            compInt.states = value
            return nothing
        end
        if key.variable == "#init_states"
            compInt.init_states = value
            return nothing
        end
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

function variables(component::DuplicatedComponent)
    return variables(component.component)
end
