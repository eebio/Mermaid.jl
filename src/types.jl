using OrdinaryDiffEqCore: ODEIntegrator

struct ConnectedVariable
    component::String
    variable::String
    variableindex::Union{Int,Nothing,UnitRange{Int}}
    fullname::String
end

function ConnectedVariable(name::AbstractString)
    return parsevariable(name)
end

struct Connector
    inputs::Vector{ConnectedVariable}
    outputs::Vector{ConnectedVariable}
    func::Union{Nothing,Function}
end

function Connector(;inputs::Vector{T}, outputs::Vector{S}, func = nothing) where T <: AbstractString where S <: AbstractString
    inputs = [ConnectedVariable(i) for i in inputs]
    outputs = [ConnectedVariable(o) for o in outputs]
    return Connector(inputs, outputs, func)
end

abstract type AbstractComponent end
# Required fields:
# - model: The model of the component
# - name: The name of the component
# - outputs: The outputs of the component
# - inputs: The inputs of the component
# - state: The state of the component

abstract type AbstractTimeIndependentComponent <: AbstractComponent end
# Required fields:
# As above, plus:
# - model: should be a function which takes inputs and state and returns outputs and edits state in place
# - precompute: Bool indicating whether to compute the solution prior to everything else

abstract type AbstractTimeDependentComponent <: AbstractComponent end

abstract type AbstractMermaidSolver end

abstract type ComponentIntegrator end

mutable struct MermaidIntegrator
    integrators::Vector
    connectors::Vector{Connector}
    maxt::Float64
    currtime::Float64
    alg::AbstractMermaidSolver
    save_vars::Vector{String}
end

@kwdef struct MermaidProblem
    components::Vector
    connectors::Vector{Connector}
    max_t::Float64 = 1.0
end

struct MermaidSolution
    t::Vector
    u::Dict
end

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

function Base.getindex(sol::MermaidSolution, var::AbstractString)
    return sol.u[parsevariable(var)]
end
