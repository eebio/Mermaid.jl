module DiffEqExt

using Mermaid
using CommonSolve
using DiffEqBase
using OrderedCollections: OrderedDict

"""
    DEComponent(model::DiffEqBase.AbstractDEProblem, alg::AbstractDEAlgorithm; name::String="DE Component",
                     state_names::Dict{String,Any}=Dict{String,Any}(),
                     time_step::Float64=1.0,
                     intkwargs::Tuple=())

# Arguments
- `model`: SciML DE problem (e.g., ODEProblem, SDEProblem, DAEProblem, etc.)
- `alg`: Algorithm from DifferentialEquations.jl to be used for solving the DEProblem.

# Keyword Arguments
- `name`: Name of the component (default: "DE Component")
- `state_names`: Dictionary mapping variable names (as strings) to their corresponding indices in the state vector or symbols from ModelingToolkit/Symbolics (default: empty dictionary)
- `time_step`: Time step for the component (default: 1.0)
- `intkwargs`: Additional keyword arguments for the DE solver (default: empty tuple)
"""
function Mermaid.DEComponent(model::DiffEqBase.AbstractDEProblem, alg::DiffEqBase.AbstractDEAlgorithm; name="DE Component",
    state_names=Dict{String,Any}(),
    time_step::Real=1.0,
    intkwargs=())
    return Mermaid.DEComponent(model, name, state_names, time_step, alg, intkwargs)
end

function CommonSolve.init(c::Mermaid.DEComponent, conns::Vector{Connector})
    integrator = Mermaid.DEComponentIntegrator(init(c.model, c.alg; dt=c.time_step, c.intkwargs...), c, OrderedDict{ConnectedVariable,Any}(), OrderedDict{ConnectedVariable,Any}())
    inputs, outputs = inputsandoutputs(integrator, conns)
    integrator.inputs = inputs
    integrator.outputs = outputs
    return integrator
end

function CommonSolve.step!(compInt::Mermaid.DEComponentIntegrator)
    for (key, value) in compInt.inputs
        setstate!(compInt, key, value)
    end
    CommonSolve.step!(compInt.integrator, compInt.component.time_step, true)
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
    return keys(component.state_names)
end

end
