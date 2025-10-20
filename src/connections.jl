"""
    ConnectedVariable <: AbstractConnectedVariable

Points to a variable within a component.

# Fields
- `component::String`: Name of the component.
- `variable::String`: Name of the variable.
- `variableindex::Union{Nothing,AbstractVector{Int},Int}`: Index or range for the variable,
    if applicable.
- `duplicatedindex::Union{Nothing,AbstractVector{Int},Int}`: Index for duplicated
    components, if applicable.
"""
struct ConnectedVariable <: AbstractConnectedVariable
    component::String
    variable::String
    variableindex::Union{Nothing, AbstractVector{Int}, Int}
    duplicatedindex::Union{Nothing, AbstractVector{Int}, Int}
end

"""
    ConnectedVariable(name::AbstractString)

Construct a [ConnectedVariable](@ref) from a string name.

# Arguments
- `name::AbstractString`: The full variable name.
"""
function ConnectedVariable(name::AbstractString)
    # Parse the variable name to extract its parts
    component, variable = split(name, ".")
    # Is there a variable index
    if contains(variable, "[")
        variable, index = split(variable, "[")
        # Strip the final "]"
        index = strip(index, ']')
        # This will parse "1:5" to a UnitRange{Int} or "3" to an Int
        index = eval(Meta.parse(index))
    else
        # No index
        index = nothing
    end
    # Is there a duplicated index
    if contains(component, "[")
        component, dupindex = split(component, "[")
        # Strip the final "]"
        dupindex = strip(dupindex, ']')
        # This will parse "1:5" to a UnitRange{Int} or "3" to an Int
        dupindex = eval(Meta.parse(dupindex))
    else
        dupindex = nothing
    end
    return ConnectedVariable(component, variable, index, dupindex)
end

"""
    Connector <: AbstractConnector

Represents a connection between multiple [ConnectedVariables](@ref ConnectedVariable),
possibly with a transformation function.

# Fields
- `inputs::Vector{<:AbstractConnectedVariable}`: Input variables for the connector.
- `outputs::Vector{<:AbstractConnectedVariable}`: Output variables for the connector.
- `func::Union{Nothing,Function}`: Optional function to transform inputs to outputs.
"""
struct Connector <: AbstractConnector
    inputs::Vector{T} where {T <: AbstractConnectedVariable}
    outputs::Vector{U} where {U <: AbstractConnectedVariable}
    func::Union{Nothing, Function}
end

"""
    Connector(inputs, outputs; func=nothing)

Construct a [Connector](@ref) from string names for inputs and outputs.

# Arguments
- `inputs::Vector{<:AbstractString}`: Names of input variables.
- `outputs::Vector{<:AbstractString}`: Names of output variables.

# Keyword Arguments
- `func`: Function for mapping inputs to outputs. Defaults to `nothing` which passes a
    single input, to all outputs.
"""
function Connector(; inputs::Vector{T}, outputs::Vector{S},
        func = nothing) where {T <: AbstractString} where {S <: AbstractString}
    inputs = [ConnectedVariable(i) for i in inputs]
    outputs = [ConnectedVariable(o) for o in outputs]
    return Connector(inputs, outputs, func)
end

"""
    fullname(var::AbstractConnectedVariable)

Return the full name of a [ConnectedVariable](@ref) as a string.

# Arguments
- `var::AbstractConnectedVariable`: The connected variable to get the full name for.

# Returns
- `String`: The full name of the connected variable.
"""
function Base.fullname(var::AbstractConnectedVariable)
    # Construct the full name of the variable
    comp = var.component
    dupindex = isnothing(var.duplicatedindex) ? "" : "[" * string(var.duplicatedindex) * "]"
    variable = var.variable
    index = isnothing(var.variableindex) ? "" : "[" * string(var.variableindex) * "]"
    return comp * dupindex * "." * variable * index
end

"""
    update_outputs!(compInt::AbstractComponentIntegrator)

Update the outputs of an [AbstractComponentIntegrator](@ref) based on its current state.

# Arguments
- `compInt::AbstractComponentIntegrator`: The component integrator whose outputs are to be
    updated.
"""
function update_outputs!(compInt::AbstractComponentIntegrator)
    # Update the outputs of the component based on the current state
    for output_key in keys(compInt.outputs)
        compInt.outputs[output_key] = getstate(compInt, output_key)
    end
end

"""
    update_inputs!(merInt::AbstractMermaidIntegrator)

Updates the inputs of each [AbstractComponentIntegrator](@ref) within the given
    [MermaidIntegrator](@ref) instance based on the outputs of other components and the
    defined [Connectors](@ref Connector).

# Arguments
- `merInt::AbstractMermaidIntegrator`: The integrator containing components, connectors,
    and their current states.
"""
function update_inputs!(merInt::AbstractMermaidIntegrator)
    # Update the inputs of the ODE component based on the outputs of other components
    for conn in merInt.connectors
        # Get the values of the connectors inputs
        inputs = []
        for input in conn.inputs
            # Find the corresponding integrator
            index = findfirst(
                i -> i.component.name == input.component, merInt.integrators)
            if index !== nothing
                integrator = merInt.integrators[index]
                # Get the value of the input from the integrator
                push!(inputs, integrator.outputs[input])
            end
        end
        if isnothing(conn.func)
            if length(inputs) == 1
                outputs = inputs[1]
            else
                outputs = inputs
            end
        else
            outputs = conn.func(inputs...)
        end
        # Set the inputs of the corresponding integrators
        for output in conn.outputs
            # Find the corresponding integrator
            index = findfirst(
                i -> i.component.name == output.component, merInt.integrators)
            if index !== nothing
                integrator = merInt.integrators[index]
                # Set the input value for the integrator
                integrator.inputs[output] = outputs
            end
        end
    end
end

"""
    inputsandoutputs(integrator::ComponentIntegrator, conns::Vector{Connector}, compName::AbstractString)

Generates the input and output dictionaries of a component integrator based on its
    connections.

# Arguments
- `integrator::AbstractComponentIntegrator`: The component integrator whose inputs and
    outputs are to be extracted.
- `conns::Vector{<:AbstractConnector}`: The connectors that define the inputs and outputs of
    the component. Extra connections not involving the `integrator` will be ignored.

# Returns
- `outputs::OrderedDict{AbstractConnectedVariable,Any}`: An ordered dictionary mapping
    [ConnectedVariable](@ref) names to their initial values from the component.
- `inputs::OrderedDict{AbstractConnectedVariable,Any}`: An ordered dictionary mapping
    [ConnectedVariable](@ref) names to their current values (initially 0).
"""
function inputsandoutputs(integrator::AbstractComponentIntegrator,
        conns::Vector{T}) where {T <: Mermaid.AbstractConnector}
    outputs = OrderedDict{AbstractConnectedVariable, Any}()
    inputs = OrderedDict{AbstractConnectedVariable, Any}()
    for conn in conns
        # If connection has an input from this component, its an output of the component
        for input in conn.inputs
            if input.component == integrator.component.name
                outputs[input] = getstate(integrator, input)
            end
        end
        # If connection has an output to this component, its an input of the component
        for output in conn.outputs
            if output.component == integrator.component.name
                inputs[output] = isnothing(output.variableindex) ? 0 :
                                 [0 for _ in output.variableindex]
            end
        end
    end
    return inputs, outputs
end
