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
