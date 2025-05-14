using CommonSolve

# Predefined concrete types
@kwdef struct PDEComponent <: AbstractTimeDependentComponent
    model::ODEProblem
    name::String = "PDE Component"
    state_names::Dict{String,Any} = Dict{String,Any}() # Dictionary that maps state names (given as strings) to their corresponding indices in the state vector (or symbols for MTK)
    time_step::Float64 = 1.0
    alg = Rodas5()
    intkwargs::Tuple{Pair{Symbol,Any}} = ()
end

mutable struct PDEComponentIntegrator <: ComponentIntegrator
    integrator::OrdinaryDiffEqCore.ODEIntegrator
    component::PDEComponent
    outputs::Dict{String,Any}
    inputs::Dict{String,Any}
end

"""
    init(c::PDEComponent)
"""
function CommonSolve.init(c::PDEComponent, conns::Vector{Connector})
    outputs = Dict{String,Any}() # Full variable name => Initial value from component
    inputs = Dict{String,Any}() # Full variable name => Value (initially 0)
    for conn in conns
        # If connection has an input from this component, store its index and function as a ComponentIntegrator.output
        for input in conn.inputs
            if split(input, ".")[1] != c.name
                continue
            end
            outputIndex = c.state_names[split(input, ".")[2]]
            # If the index is a MTK symbol then get the variable index
            if !isa(outputIndex, Integer)
                outputIndex = variable_index(c.model.f.sys, outputIndex)
            end
            outputs[input] = c.model.u0[outputIndex]
        end
        for output in conn.outputs
            if split(output, ".")[1] != c.name
                continue
            end
            inputs[output] = 0
        end
    end
    integrator = PDEComponentIntegrator(init(c.model, c.alg; dt=c.time_step, c.intkwargs...), c, outputs, inputs)
    return integrator
end

"""
    step!(compInt::PDEComponentIntegrator)

Steps the PDE component integrator.
# Arguments
- `compInt::PDEComponentIntegrator`: The PDE component integrator to be stepped.
"""
function CommonSolve.step!(compInt::PDEComponentIntegrator)
    for (key, value) in compInt.inputs
        setstate!(compInt, split(key, ".")[2], value)
    end
    u_modified!(compInt.integrator, true)
    CommonSolve.step!(compInt.integrator)
end
