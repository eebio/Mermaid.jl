module AgentsExt

using Mermaid
using CommonSolve
using Agents
using OrderedCollections: OrderedDict

export AgentsComponent

"""
    AgentsComponent(model::StandardABM; name="Agents Component",
                    state_names=Dict{String,Any}(), time_step::Real=1.0)

A Mermaid component that wraps an agent-based model (ABM) using the Agents.jl package.

# Arguments
- `model::StandardABM`: The agent-based model to be solved.

# Keyword Arguments
- `name::AbstractString`: The name of the component. Defaults to "Agents Component".
- `state_names`: A dictionary mapping variable names (as strings) to their corresponding
    properties (agent properties or model properties) in the `model`. Defaults to an empty
    dictionary.
- `time_step::Real`: The time step for the component (not the ABM solver timestep), i.e. how
    frequently should the inputs and outputs be updated.
"""
function Mermaid.AgentsComponent(model::StandardABM;
        name::AbstractString = "Agents Component", state_names = Dict{String, Any}(),
        time_step::Real = 1.0)
    return Mermaid.AgentsComponent(model, name, state_names, time_step)
end

function CommonSolve.init(c::AgentsComponent)
    integrator = AgentsComponentIntegrator(deepcopy(c.model), c)
    return integrator
end

function CommonSolve.step!(compInt::AgentsComponentIntegrator)
    step!(compInt.integrator)
end

function Mermaid.getstate(compInt::AgentsComponentIntegrator, key::ConnectedVariable)
    if first(key.variable) == '#'
        # Special variables
        if key.variable == "#time"
            return abmtime(compInt.integrator) * compInt.component.time_step
        end
        if key.variable == "#model"
            return getstate(compInt; copy = true)
        end
        if key.variable == "#ids"
            return collect(allids(compInt.integrator))
        end
    end
    index = compInt.component.state_names[key.variable]
    if isnothing(key.variableindex)
        # If model level property exists, return it directly
        !isnothing(abmproperties(compInt.integrator)) &&
            haskey(abmproperties(compInt.integrator), index) &&
            return getproperty(compInt.integrator, index)
        # Otherwise, assume it's an agent property and return it for all agents
        return [getproperty(i, index) for i in allagents(compInt.integrator)]
    else
        # If model level property exists, return it after indexing
        !isnothing(abmproperties(compInt.integrator)) &&
            haskey(abmproperties(compInt.integrator), index) &&
            return getproperty(compInt.integrator, index)[key.variableindex]
        # Otherwise, assume it's an agent property and return it for all agents in range
        out = [getproperty(compInt.integrator[i], index) for i in key.variableindex]
        if length(out) == 1
            return out[1] # If only one agent, return the single value
        else
            return out # Otherwise, return the vector of values
        end
    end
end

function Mermaid.setstate!(
        compInt::AgentsComponentIntegrator, key::ConnectedVariable, value)
    # TODO add ids exception? Is there a way to specify the id when creating an agent
    if first(key.variable) == '#'
        # Special variables
        if key.variable == "#time"
            @warn "Agents.jl does not support setting time directly. The time is stored within #model. DuplicatedComponents of Agent-based models still works."
            return nothing
        end
        if key.variable == "#model"
            setstate!(compInt, value)
            return nothing
        end
    end
    index = compInt.component.state_names[key.variable]
    if isnothing(key.variableindex)
        # If model level property exists, return it directly
        if !isnothing(abmproperties(compInt.integrator)) &&
           haskey(abmproperties(compInt.integrator), index)
            setindex!(abmproperties(compInt.integrator), value, index)
        else
            # Otherwise, assume it's an agent property and set it for all agents
            for (i, agent) in enumerate(allagents(compInt.integrator))
                setproperty!(agent, index, value[i])
            end
        end
    else
        # If model level property exists, return it after indexing
        if !isnothing(abmproperties(compInt.integrator)) &&
           haskey(abmproperties(compInt.integrator), index)
            k = 1
            for i in key.variableindex
                abmproperties(compInt.integrator)[index][i] = value[k]
                k += 1
            end
        else
            # Otherwise, assume it's an agent property and set it for all agents in range
            k = 1
            for i in key.variableindex
                setproperty!(compInt.integrator[i], index, value[k])
                k += 1
            end
        end
    end
    return nothing
end

function Mermaid.getstate(compInt::AgentsComponentIntegrator)
    return compInt.integrator
end

function Mermaid.setstate!(compInt::AgentsComponentIntegrator, state::StandardABM)
    compInt.integrator = state
end

function Mermaid.variables(component::AgentsComponent)
    return union(keys(component.state_names), ["#model", "#time"])
end
end
