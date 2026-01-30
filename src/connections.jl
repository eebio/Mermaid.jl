"""
    ConnectedVariable <: AbstractConnectedVariable

Points to a variable within a component.

# Fields
- `component::String`: Name of the component.
- `variable::String`: Name of the variable.
- `variableindex::Union{Nothing,Vector{Int},Int}`: Index or range for the variable,
    if applicable.
- `duplicatedindex::Union{Nothing,Vector{Int},Int}`: Index for duplicated
    components, if applicable.
"""
struct ConnectedVariable{U <: AbstractString, V <: AbstractString, X <: Union{Nothing, AbstractVector{Int}, Int}, Y <: Union{Nothing, AbstractVector{Int}, Int}} <: AbstractConnectedVariable
    component::U
    variable::V
    variableindex::X
    duplicatedindex::Y
end

function Base.hash(cv::ConnectedVariable, h::UInt) return hash(fullname(cv), h) end
function Base.isequal(cv1::ConnectedVariable, cv2::ConnectedVariable)
    return fullname(cv1) == fullname(cv2)
end

"""
    ConnectedVariable(name::AbstractString)

Construct a [ConnectedVariable](@ref) from a string name.

# Arguments
- `name::AbstractString`: The full variable name.

# Examples
`ConnectedVariable("comp.var")`
`ConnectedVariable("comp.var[1:5]")` where 1:5 is the variableindex
`ConnectedVariable("comp[2].var")` where 2 is the duplicatedindex
`ConnectedVariable("comp[1:3].var[4]")`
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
    inputs::Vector{ConnectedVariable}
    outputs::Vector{ConnectedVariable}
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
    runconnection(merInt::AbstractMermaidIntegrator, conn::AbstractConnector)

Extract all the input states from `merInt`, apply the connection function, and return the
output.

# Arguments
- `merInt::AbstractMermaidIntegrator`: The Mermaid integrator containing the components.
- `conn::AbstractConnector`: The connector defining the connection.
"""
function runconnection(merInt::AbstractMermaidIntegrator, conn::AbstractConnector)
    # Get the values of the connectors inputs
    inputs = []
    for input in conn.inputs
        # Find the corresponding integrator
        index = findfirst(
            i -> name(i) == input.component, merInt.integrators)
        if index !== nothing
            integrator = merInt.integrators[index]
            # Get the value of the input from the integrator
            push!(inputs, getstate(integrator, input))
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
    return outputs
end

"""
    runconnection!(merInt::AbstractMermaidIntegrator, conn::AbstractConnector)

Extract all the input states from `merInt`, apply the connection function, and set the output
states in `merInt`.

# Arguments
- `merInt::AbstractMermaidIntegrator`: The Mermaid integrator containing the components.
- `conn::AbstractConnector`: The connector defining the connection.
"""
function runconnection!(merInt::AbstractMermaidIntegrator, conn::AbstractConnector)
    outputs = runconnection(merInt, conn)
    # Set the outputs in the corresponding integrators
    for output in conn.outputs
        # Find the corresponding integrator
        index = findfirst(
            i -> name(i) == output.component, merInt.integrators)
        if index !== nothing
            integrator = merInt.integrators[index]
            # Set the input value for the integrator
            setstate!(integrator, output, outputs)
        end
    end
end
