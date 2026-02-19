module DiffEqExt

using Mermaid
using CommonSolve
using DiffEqBase
using OrderedCollections: OrderedDict

"""
    DEComponent(model::DiffEqBase.AbstractDEProblem, alg;
                name::String="DE Component", time_step::Float64=1.0, intkwargs::Tuple=(),
                state_names::Dict{String,Any}=Dict{String,Any}())
    DEComponent(model::DiffEqBase.AbstractDEProblem; kwargs...)

# Arguments
- `model::DiffEqBase.AbstractDEProblem`: SciML DE problem (e.g., ODEProblem, etc.)
- `alg`: Algorithm from DifferentialEquations.jl to be used
    for solving the DEProblem. If no algorithm is provided, the algorithm will be
    automatically chosen by DifferentialEquations.jl.

# Keyword Arguments
- `name::AbstractString`: Name of the component. Defaults to "DE Component".
- `time_step::Real`: Time step for the component. Defaults to 1.0.
- `intkwargs`: Additional keyword arguments for the DE solver. Defaults to no keywords.
- `state_names`: Dictionary mapping variable names (as strings) to their corresponding
    indices in the state vector or symbols from Symbolics.jl. Defaults to an empty
    dictionary.
"""
function Mermaid.DEComponent(model::DiffEqBase.AbstractDEProblem,
        alg; name = "DE Component", time_step::Real = 1.0, intkwargs = (),
        state_names = Dict{String, Any}())
    return Mermaid.DEComponent(model, name, state_names, time_step, alg, intkwargs)
end

function Mermaid.DEComponent(model::DiffEqBase.AbstractDEProblem; kwargs...)
    return Mermaid.DEComponent(model, nothing; kwargs...)
end

function CommonSolve.init(c::Mermaid.DEComponent)
    integrator = Mermaid.DEComponentIntegrator(
        init(c.model, c.alg; c.intkwargs...), c)
    return integrator
end

function CommonSolve.step!(compInt::Mermaid.DEComponentIntegrator)
    CommonSolve.step!(compInt.integrator, time_step(compInt), true)
end

function Mermaid.getstate(compInt::Mermaid.DEComponentIntegrator, key)
    if first(key.variable) == '#'
        if key.variable == "#time"
            return compInt.integrator.t
        end
    end
    index = compInt.component.state_names[key.variable]
    return compInt.integrator[index]
end

function Mermaid.getstate(compInt::Mermaid.DEComponentIntegrator)
    return compInt.integrator.u
end

function Mermaid.setstate!(compInt::Mermaid.DEComponentIntegrator, key, value)
    u_modified!(compInt.integrator, true)
    if first(key.variable) == '#'
        if key.variable == "#time"
            u_modified!(compInt.integrator, true)
            compInt.integrator.t = value
            return nothing
        end
    end
    index = compInt.component.state_names[key.variable]
    compInt.integrator[index] = value
end

function Mermaid.setstate!(compInt::Mermaid.DEComponentIntegrator, value)
    u_modified!(compInt.integrator, true)
    compInt.integrator.u = value
end

function Mermaid.variables(component::Mermaid.DEComponent)
    return union(keys(component.state_names), ["#time"])
end

end
