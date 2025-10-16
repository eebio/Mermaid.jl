module MethodOfLinesExt

using Mermaid
using CommonSolve
using OrderedCollections: OrderedDict
using DiffEqBase

"""
    MOLComponent(model::DiffEqBase.AbstractDEProblem, alg::DiffEqBase.AbstractDEAlgorithm;
                 name::String="MOL Component", time_step::Real=1.0, intkwargs::Tuple=(),
                 state_names::Dict{String,Any}=Dict{String,Any}())

# Arguments
- `model::DiffEqBase.AbstractDEProblem`: SciML DE problem (e.g., ODEProblem, etc.)
- `alg::DiffEqBase.AbstractDEAlgorithm`: Algorithm from DifferentialEquations.jl to be used
    for solving the DEProblem.

# Keyword Arguments
- `name::AbstractString`: Name of the component. Defaults to "MOL Component".
- `time_step::Real`: Time step for the component. Defaults to 1.0.
- `intkwargs`: Additional keyword arguments for the DE solver. Defaults to no keywords.
- `state_names`: Dictionary mapping variable names (as strings) to their corresponding
    indices in the state vector or symbols from Symbolics.jl. Defaults to an empty
    dictionary.
"""
function Mermaid.MOLComponent(model::DiffEqBase.AbstractDEProblem,
        alg::DiffEqBase.AbstractDEAlgorithm; name = "MOL Component",
        time_step::Real = 1.0, intkwargs = (), state_names = Dict{String, Any}())
    return Mermaid.MOLComponent(model, name, state_names, time_step, alg, intkwargs)
end

function CommonSolve.init(
        c::MOLComponent, conns::Vector{T}) where {T <: Mermaid.AbstractConnector}
    integrator = MOLComponentIntegrator(
        init(c.model, c.alg; dt = c.time_step, c.intkwargs...), c,
        OrderedDict{ConnectedVariable, Any}(), OrderedDict{ConnectedVariable, Any}())
    inputs, outputs = inputsandoutputs(integrator, conns)
    integrator.inputs = inputs
    integrator.outputs = outputs
    return integrator
end

function CommonSolve.step!(compInt::MOLComponentIntegrator)
    for (key, value) in compInt.inputs
        setstate!(compInt, key, value)
    end
    CommonSolve.step!(compInt.integrator, compInt.component.time_step, true)
end

function Mermaid.getstate(compInt::MOLComponentIntegrator, key)
    if first(key.variable) == '#'
        if key.variable == "#time"
            return compInt.integrator.t
        end
    end
    if isnothing(key.variableindex)
        # No index for variable
        index = compInt.component.state_names[key.variable]
        return compInt.integrator[index]
    else
        index = compInt.component.state_names[key.variable]
        return compInt.integrator[index][key.variableindex]
    end
end

function Mermaid.getstate(compInt::MOLComponentIntegrator)
    # Return the full state vector
    return compInt.integrator.u
end

function Mermaid.setstate!(compInt::MOLComponentIntegrator, key, value)
    u_modified!(compInt.integrator, true)
    if first(key.variable) == '#'
        if key.variable == "#time"
            compInt.integrator.t = value
            return nothing
        end
    end
    if isnothing(key.variableindex)
        # No index for variable
        index = compInt.component.state_names[key.variable]
        compInt.integrator[index] = value
    else
        key.variableindex isa AbstractVector
        # Single index for one variable
        index = compInt.component.state_names[key.variable]
        if length(index[key.variableindex]) == 1
            # If the index is a single value, set it directly
            compInt.integrator.u[(index)[key.variableindex][1]] = value
        else
            # If the index is a vector, set the corresponding values
            compInt.integrator.u[(index)[key.variableindex]] .= value
        end
    end
end

function Mermaid.setstate!(compInt::MOLComponentIntegrator, value)
    # Set the full state vector
    compInt.integrator.u = value
end

function Mermaid.variables(component::MOLComponent)
    return union(keys(component.state_names), ["#time"])
end

end
