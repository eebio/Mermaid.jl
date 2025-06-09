using OrdinaryDiffEqCore: ODEIntegrator

"""
    ConnectedVariable

Points to a variable that is connected between components.

# Fields
- `component::String`: Name of the component.
- `variable::String`: Name of the variable.
- `variableindex::Union{Nothing,AbstractVector{Int},Int}`: Index or range for the variable, if applicable.
- `duplicatedindex::Union{Nothing,AbstractVector{Int},Int}`: Index for duplicated components, if applicable.
"""
struct ConnectedVariable
    component::String
    variable::String
    variableindex::Union{Nothing,AbstractVector{Int},Int}
    duplicatedindex::Union{Nothing,AbstractVector{Int},Int}
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
struct Connector
    inputs::Vector{ConnectedVariable}
    outputs::Vector{ConnectedVariable}
    func::Union{Nothing,Function}
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
function Connector(; inputs::Vector{T}, outputs::Vector{S}, func=nothing) where T<:AbstractString where S<:AbstractString
    inputs = [ConnectedVariable(i) for i in inputs]
    outputs = [ConnectedVariable(o) for o in outputs]
    return Connector(inputs, outputs, func)
end

abstract type AbstractComponent end

abstract type AbstractTimeIndependentComponent <: AbstractComponent end

abstract type AbstractTimeDependentComponent <: AbstractComponent end

"""
    AbstractMermaidSolver
Base type for all solvers in the Mermaid framework.
"""
abstract type AbstractMermaidSolver end

"""
    ComponentIntegrator
Base type for all component integrators in the Mermaid framework.
"""
abstract type ComponentIntegrator end

"""
    MermaidIntegrator

Mutable struct for integrating a hybrid [MermaidProblem](@ref).
This struct holds all the [ComponentIntegrators](@ref ComponentIntegrator) and [Connectors](@ref Connector) to store the current state of the hybrid simulation.

# Fields
- `integrators::Vector`: Vector of [ComponentIntegrator](@ref).
- `connectors::Vector{Connector}`: Vector of [Connector](@ref).
- `maxt::Float64`: Maximum simulation time.
- `currtime::Float64`: Current simulation time.
- `alg::AbstractMermaidSolver`: Algorithm used for integration.
- `save_vars::Vector{String}`: Variables to save during integration.
"""
mutable struct MermaidIntegrator
    integrators::Vector
    connectors::Vector{Connector}
    maxt::Float64
    currtime::Float64
    alg::AbstractMermaidSolver
    save_vars::Vector{String}
end

"""
    MermaidProblem

Struct for defining a Mermaid problem.
This struct contains the components, connectors and other properties of the hybrid simulation.

# Fields
- `components::Vector`: Vector of Components.
- `connectors::Vector{Connector}`: Vector of [Connector](@ref).
- `max_t::Float64`=1.0: Maximum simulation time.
"""
@kwdef struct MermaidProblem
    components::Vector
    connectors::Vector{Connector}
    max_t::Float64 = 1.0
end

"""
    MermaidSolution

Struct for storing the solution of a hybrid Mermaid simulation.

# Fields
- `t::Vector`: Time points.
- `u::Dict`: Dictionary mapping variables to their solution arrays.
"""
struct MermaidSolution
    t::Vector
    u::Dict
end

"""
    MermaidSolution(int::MermaidIntegrator)

Construct a [MermaidSolution](@ref) from a [MermaidIntegrator](@ref).

# Arguments
- `int::MermaidIntegrator`: The integrator to extract solution structure from.

# Returns
- `MermaidSolution`: The initialized solution object.
"""
function MermaidSolution(int::MermaidIntegrator)
    u = Dict()
    # TODO: This is still lacking, if int.save_vars is comp.u[5] but state_names only describes comp.u, this won't work
    for i in int.integrators
        for key in keys(i.component.state_names)
            fullname = join([i.component.name, key], ".")
            if length(int.save_vars) == 0 || fullname in int.save_vars
                u[parsevariable(fullname)] = []
            end
        end
    end
    return MermaidSolution([], u)
end

"""
    Base.getindex(sol::MermaidSolution, var::AbstractString)

Get the solution array for a variable from a [MermaidSolution](@ref).

# Arguments
- `sol::MermaidSolution`: The solution object.
- `var::AbstractString`: The variable name.

# Returns
- The solution array for the specified variable.
"""
function Base.getindex(sol::MermaidSolution, var::AbstractString)
    var = parsevariable(var)
    if haskey(sol.u, var)
        return sol.u[var]
    else
        # See if we have a key without an index
        var_no_index = ConnectedVariable(var.component, var.variable, nothing, nothing) # TODO I'm not sure how the duplicatedindex data is stored in the solution
        return [i[var.variableindex] for i in sol.u[var_no_index]]
    end
end
