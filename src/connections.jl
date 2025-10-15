"""
    ConnectedVariable

Points to a variable that is connected between components.

# Fields
- `component::String`: Name of the component.
- `variable::String`: Name of the variable.
- `variableindex::Union{Nothing,AbstractVector{Int},Int}`: Index or range for the variable, if applicable.
- `duplicatedindex::Union{Nothing,AbstractVector{Int},Int}`: Index for duplicated components, if applicable.
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

# Returns
- `ConnectedVariable`: The parsed connected variable.
"""
function ConnectedVariable(name::AbstractString)
    return parsevariable(name)
end

"""
    Connector

Represents a connection between multiple [ConnectedVariables](@ref ConnectedVariable), possibly with a transformation function.

# Fields
- `inputs::Vector{ConnectedVariable}`: Input variables for the connector.
- `outputs::Vector{ConnectedVariable}`: Output variables for the connector.
- `func::Union{Nothing,Function}`: Optional function to transform inputs to outputs.
"""
struct Connector <: AbstractConnector
    inputs::Vector{AbstractConnectedVariable}
    outputs::Vector{AbstractConnectedVariable}
    func::Union{Nothing, Function}
end

"""
    Connector(;inputs, outputs, func=nothing)

Construct a [Connector](@ref) from string names for inputs and outputs.

# Keyword Arguments
- `inputs::Vector{<:AbstractString}`: Names of input variables.
- `outputs::Vector{<:AbstractString}`: Names of output variables.
- `func`: Optional function for transformation (default: `nothing`).

# Returns
- `Connector`: The constructed connector.
"""
function Connector(; inputs::Vector{T}, outputs::Vector{S},
        func = nothing) where {T <: AbstractString} where {S <: AbstractString}
    inputs = [ConnectedVariable(i) for i in inputs]
    outputs = [ConnectedVariable(o) for o in outputs]
    return Connector(inputs, outputs, func)
end

"""
    parsevariable(name::AbstractString) -> ConnectedVariable

Parses a variable name as a foramtted string to a [ConnectedVariable](@ref).

# Arguments
- `name::AbstractString`: The variable name to parse. It can include an optional index, which may be a single integer (e.g., `"comp.var[3]"`) or a range (e.g., `"comp.var[1:5]"`).

# Returns
- `ConnectedVariable`: The corresponding [ConnectedVariable](@ref), containing a `component`, `variable`, `index` (which can be `nothing`, an `Int`, or a `UnitRange{Int}`), and the original `name`.
"""
function parsevariable(name)
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

function Base.fullname(var::AbstractConnectedVariable)
    # Construct the full name of the variable
    comp = var.component
    dupindex = isnothing(var.duplicatedindex) ? "" : "[" * string(var.duplicatedindex) * "]"
    variable = var.variable
    index = isnothing(var.variableindex) ? "" : "[" * string(var.variableindex) * "]"
    return comp * dupindex * "." * variable * index
end

"""
    update_outputs!(compInt::ComponentIntegrator)

Update the outputs field of a [ComponentIntegrator](@ref) based on its current state.

# Arguments
- `compInt::ComponentIntegrator`: The component integrator whose outputs are to be updated.
"""
function update_outputs!(compInt::AbstractComponentIntegrator)
    # Update the outputs of the component based on the current state
    for output_key in keys(compInt.outputs)
        compInt.outputs[output_key] = getstate(compInt, output_key)
    end
end

"""
    update_inputs!(mermaidInt::MermaidIntegrator)

Updates the input values of each [ComponentIntegrator](@ref) within the given [MermaidIntegrator](@ref) instance based on the outputs of other components and the defined [Connectors](@ref Connector).

# Arguments
- `mermaidInt::MermaidIntegrator`: The integrator containing components, connectors, and their current states.

# Description
For each connector in `mermaidInt.connectors`, this function:
- Collects the outputs from the source components specified in the connector's inputs.
- If a function (`conn.func`) is defined for the connector, applies it to the collected inputs to compute the outputs; otherwise, passes the inputs directly.
- Assigns the resulting outputs to the appropriate input fields of the destination components specified in the connector's outputs.
"""
function update_inputs!(mermaidInt::AbstractMermaidIntegrator)
    # Update the inputs of the ODE component based on the outputs of other components
    for conn in mermaidInt.connectors
        # Get the values of the connectors inputs
        inputs = []
        for input in conn.inputs
            # Find the corresponding integrator
            index = findfirst(
                i -> i.component.name == input.component, mermaidInt.integrators)
            if index !== nothing
                integrator = mermaidInt.integrators[index]
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
                i -> i.component.name == output.component, mermaidInt.integrators)
            if index !== nothing
                integrator = mermaidInt.integrators[index]
                # Set the input value for the integrator
                integrator.inputs[output] = outputs
            end
        end
    end
end

"""
    inputsandoutputs(integrator::ComponentIntegrator, conns::Vector{Connector}, compName::AbstractString)

Generates the inputs and outputs of a component integrator based on its connections.

# Arguments
- `integrator::ComponentIntegrator`: The component integrator whose inputs and outputs are to be extracted.
- `conns::Vector{Connector}`: The connectors that define the inputs and outputs of the component.

# Returns
- `outputs::OrderedDict{ConnectedVariable,Any}`: An ordered dictionary mapping [ConnectedVariable](@ref) names to their initial values from the component.
- `inputs::OrderedDict{ConnectedVariable,Any}`: An ordered dictionary mapping [ConnectedVariable](@ref) names to their current values (initially 0).
"""
function inputsandoutputs(integrator::AbstractComponentIntegrator,
        conns::Vector{T}) where {T <: Mermaid.AbstractConnector}
    outputs = OrderedDict{AbstractConnectedVariable, Any}() # Full variable name => Initial value from component
    inputs = OrderedDict{AbstractConnectedVariable, Any}() # Full variable name => Value (initially 0)
    for conn in conns
        # If connection has an input from this component, store its index and function as a ComponentIntegrator.output
        for input in conn.inputs
            if input.component == integrator.component.name
                outputs[input] = getstate(integrator, input)
            end
        end
        for output in conn.outputs
            if output.component == integrator.component.name
                inputs[output] = isnothing(output.variableindex) ? 0 :
                                 [0 for _ in output.variableindex]
            end
        end
    end
    return inputs, outputs
end
