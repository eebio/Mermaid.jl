using CommonSolve

@kwdef struct PDEComponent <: AbstractTimeDependentComponent
    model::ODEProblem
    name::String = "PDE Component"
    state_names::Dict{String,Any} = Dict{String,Any}() # Dictionary that maps state names (given as strings) to their corresponding indices in the state vector (or symbols for MTK)
    time_step::Float64 = 1.0
    alg = Rodas5()
    intkwargs::Tuple = ()
end

mutable struct PDEComponentIntegrator <: ComponentIntegrator
    integrator::OrdinaryDiffEqCore.ODEIntegrator
    component::PDEComponent
    outputs::Dict{ConnectedVariable,Any}
    inputs::Dict{ConnectedVariable,Any}
end

function CommonSolve.init(c::PDEComponent, conns::Vector{Connector})
    outputs = Dict{ConnectedVariable,Any}() # Full variable name => Initial value from component
    inputs = Dict{ConnectedVariable,Any}() # Full variable name => Value (initially 0)
    for conn in conns
        # If connection has an input from this component, store its index and function as a ComponentIntegrator.output
        for input in conn.inputs
            if input.component == c.name
                outputIndex = c.state_names[input.variable]
                # If the index is a MTK symbol then get the variable index
                if isa(outputIndex, Num) || isa(outputIndex, Symbolics.Arr)
                    outputIndex = variable_index(c.model.f.sys, outputIndex)
                end
                # TODO I think we can allow variable indexes here, for if each element of the state is a vector
                outputs[input] = c.model.u0[outputIndex]
            end
        end
        for output in conn.outputs
            if output.component == c.name
                inputs[output] = 0
            end
        end
    end
    integrator = PDEComponentIntegrator(init(c.model, c.alg; dt=c.time_step, c.intkwargs...), c, outputs, inputs)
    return integrator
end

function CommonSolve.step!(compInt::PDEComponentIntegrator)
    for (key, value) in compInt.inputs
        setstate!(compInt, key, value)
    end
    u_modified!(compInt.integrator, true)
    CommonSolve.step!(compInt.integrator, compInt.component.time_step, true)
end

function getstate(compInt::PDEComponentIntegrator, key)
    if isnothing(key.variableindex)
        # No index for variable
        index = compInt.component.state_names[key.variable]
        return compInt.integrator[index]
    else
        index = compInt.component.state_names[key.variable]
        return compInt.integrator[index][key.variableindex]
    end
end

function getstate(compInt::PDEComponentIntegrator)
    # Return the full state vector
    return compInt.integrator.u
end

function setstate!(compInt::PDEComponentIntegrator, key, value)
    if isnothing(key.variableindex)
        # No index for variable
        index = compInt.component.state_names[key.variable]
        compInt.integrator[index] = value
    else key.variableindex isa AbstractVector
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

function setstate!(compInt::PDEComponentIntegrator, value)
    # Set the full state vector
    compInt.integrator.u = value
end

function gettime(compInt::PDEComponentIntegrator)
    return compInt.integrator.t
end

function settime!(compInt::PDEComponentIntegrator, t)
    u_modified!(compInt.integrator, true)
    compInt.integrator.t = t
end
