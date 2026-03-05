module DiffEqExt

using Mermaid
using CommonSolve
using DiffEqBase
using JumpProcesses
using OrderedCollections: OrderedDict

"""
    JumpComponent(model::JumpProblem, alg;
                name::String="Jump Component", timestep::Float64=1.0, intkwargs::Tuple=(),
                state_names::Dict{String,Any}=Dict{String,Any}())

# Arguments
- `model::JumpProblem`: SciML Jump problem.
- `alg`: Algorithm from DifferentialEquations.jl to be used
    for solving the JumpProblem.

# Keyword Arguments
- `name::AbstractString`: Name of the component. Defaults to "Jump Component".
- `timestep::Real`: Time step for the component. Defaults to 1.0.
- `intkwargs`: Additional keyword arguments for the Jump solver. Defaults to no keywords.
- `state_names`: Dictionary mapping variable names (as strings) to their corresponding
    indices in the state vector or symbols from Symbolics.jl. Defaults to an empty
    dictionary.
"""
function Mermaid.JumpComponent(model::JumpProblem,
        alg; name = "Jump Component", timestep::Real = 1.0, intkwargs = (),
        state_names = Dict{String, Any}())
    return Mermaid.JumpComponent(model, name, state_names, timestep, alg, intkwargs)
end

function CommonSolve.init(c::Mermaid.JumpComponent)
    integrator = Mermaid.JumpComponentIntegrator(
        init(c.model, c.alg; c.intkwargs...), c)
    return integrator
end

function CommonSolve.step!(compInt::Mermaid.JumpComponentIntegrator)
    CommonSolve.step!(compInt.integrator, timestep(compInt), true)
end

function Mermaid.getstate(compInt::Mermaid.JumpComponentIntegrator, key)
    if first(key.variable) == '#'
        if key.variable == "#time"
            return compInt.integrator.t
        end
    end
    index = compInt.component.state_names[key.variable]
    return compInt.integrator[index]
end

function Mermaid.getstate(compInt::Mermaid.JumpComponentIntegrator)
    return compInt.integrator.u
end

function Mermaid.setstate!(compInt::Mermaid.JumpComponentIntegrator, key, value)
    u_modified!(compInt.integrator, true)
    # Clear jump caches too
    # If parameters have changed, we need to do extra work
    param_changed = false
    reset_aggregated_jumps!(compInt.integrator; update_jump_params = param_changed)
    if first(key.variable) == '#'
        if key.variable == "#time"
            compInt.integrator.t = value
            return nothing
        end
    end
    index = compInt.component.state_names[key.variable]
    compInt.integrator[index] = value
end

function Mermaid.setstate!(compInt::Mermaid.JumpComponentIntegrator, value)
    u_modified!(compInt.integrator, true)
    param_changed = false
    reset_aggregated_jumps!(compInt.integrator; update_jump_params = param_changed)
    compInt.integrator.u = value
end

function Mermaid.variables(component::Mermaid.JumpComponent)
    return union(keys(component.state_names), ["#time"])
end

end
